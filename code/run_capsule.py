"""Reproducible-run entry for the automated ophys (2p) <-> Xenium coregistration capsule.

Wraps the `xenium_autocoreg` package from
https://github.com/jkim0731/2p2xenium (installed in environment/postInstall) for one
subject. The whole pipeline is GT-free: anchor-section selection, automatic pose
seeding (soma-print point-cloud descriptor), 3D tilt fitting, chain-refine propagation
across sections, fine mask-TPS registration + probability-filtered cell matching, and
3D point mapping of every Xenium cell into z-stack coordinates.

`xenium_autocoreg` itself is data-layout-agnostic (it only knows about
`xenium_autocoreg.config.SubjectConfig`) -- `subject_resolver.py` (this capsule's own sibling
module) is the CodeOcean/company-specific glue that resolves one subject's mounted assets into a
`SubjectConfig`, using this lab's asset-naming convention.

Pipeline stages (see the 2p2xenium README for full detail; all run via
``xenium_autocoreg.cli.run_subject``):
    1. anchor-section selection      2. pose seeding (auto / center-rotation)
    3. tilt fitting                  4. chain-refine propagation (all sections)
    5. fine registration + cell matching (per section)
    6. 3D point mapping of every Xenium cell centroid

Inputs (attach as data assets, mounted under ``--input_dir`` = /root/capsule/data; resolved
by this capsule's own ``subject_resolver.resolve_subject`` -- glob, no IDs hard-coded):
    Xenium-ophys-coregistered_{sid}_*      aligned-frame Xenium (Xenium_images/,
                                            Xenium_segmentation_masks/, section-to-section
                                            aligned -- see the 2p2xenium README)
    ophys-z-stacks_{sid}_registered_*      z-stack registered intensity volume (+ its
    ophys-z-stacks_{sid}_segmented_*       matching segmentation label volume), OR the
                                            multiplane-ophys_{sid}_*cortical-zstack-{registration,
                                            segmentation}_* naming (auto-detected fallback)
    Xenium_{sid}_*_processed               OPTIONAL: reporter-transcript population source
                                            (SYFP2/EGFP); falls back to all segmented Xenium
                                            cells if not attached

Outputs (to ``/root/capsule/results``, captured as a derived data asset), written directly
by ``xenium_autocoreg.cli.run_subject`` in the layout documented in the 2p2xenium README
(Affine matrices/, Xenium_affine_transformed/, warped_zstacks/, post_affine_warping/,
cell_matching_probability/, cell_centroids/, mapped_3d_coordinates/, ophys-z-stacks*/, QC/,
propagation_summary.json), plus this capsule's own ``coreg_manifest.json`` (run
parameters + per-section summary, for provenance). The package's own transient working
files (``_chain_internal/``) are moved to ``/root/capsule/scratch`` -- regenerable, not a
scientific output (see the 2p2xenium README, which labels it "harmless scratch").

Pose-seeding protocols (``--pose_mode``) -- parameter names match ``xenium_autocoreg.cli``'s own
CLI / ``--pose-json`` shape exactly (see the 2p2xenium README), so the same ``--pose_json`` file
works unchanged against either tool:
    auto             -- ``--anchor_sec`` optional (blank = auto-select).
    center-rotation  -- ``--center_um`` + ``--rotation_deg`` required.

``--zstack_scale_to_xenium`` OPTIONALLY overrides the subject's own
``SubjectConfig.zstack_scale_to_Xenium`` (the z-stack-to-Xenium linear PHYSICAL scale factor used
to seed the initial pose search) for this run -- applies to EVERY ``--pose_mode`` alike (auto and
center-rotation both ultimately seed from ``cfg.zstack_scale_to_Xenium``), not just
center-rotation. This is an INPUT prior, not the same thing as the per-candidate FITTED affine
scale reported in 2p2xenium's own logs/QC (e.g. a grid-search candidate line's ``scale=0.809``) --
that's a measured OUTPUT of the registration, not a knob.

(A third, corner-based mode is a documented TODO in ``xenium_autocoreg.pose_seed.seed_from_corners``
-- not yet implemented upstream, and not exposed here.)

``--zstack_registered_tif``/``--zstack_segmented_tif`` OPTIONALLY pin the exact z-stack files to
use, bypassing ``subject_resolver``'s automatic ``ophys-z-stacks_*``/``multiplane-ophys_*``
discovery entirely for this run. Use this when the registration/segmentation pair you want isn't
(or can't be) auto-discovered -- e.g. a differently-structured derived asset (no ``channel_0_ref_0``
nesting) that the resolver's fixed globs don't see, or you want a specific channel/sub-acquisition
the max-ROI-count rule wouldn't pick. Both must be given together (a mismatched auto+pinned pair
would silently combine two different physical acquisitions). Paths are resolved relative to
``--input_dir`` if not absolute. ``aligned_dir``/``reporter_zarr_root`` are still resolved
automatically as usual; the automatic z-stack discovery itself (``_find_zstack_pair``/
``_find_zstack_pair_multiplane``) is skipped ENTIRELY when both are given -- no
``ophys-z-stacks_*``/``multiplane-ophys_*`` asset needs to be auto-discoverable, or even attached,
for this run. Calibration (``zstack_xy_um``) is auto-derived from a ``roi_groups_metadata.json``
found near the pinned registered tif itself (not from any auto-discovered stack), falling back to
a nominal 700um/512px guess (with a printed warning) if none is found there either.
"""
import argparse
import json
import shutil
import sys
from pathlib import Path

import subject_resolver


_COORD_FLAGS = ("--center_um",)


def _fix_negative_coord_tokens(argv):
    """Rewrite `--flag value` to `--flag=value` for the comma-joined coordinate flags, so a value
    starting with a minus sign (e.g. '-50,-100') isn't misparsed by argparse as a new option --
    argparse's own negative-number heuristic only recognizes a bare `-123`/`-1.5`, not a
    comma-separated list, so `--center_um -50,-100` (exactly how CodeOcean's app panel invokes a
    capsule: `--<key> <value>` as separate argv tokens) would otherwise fail with 'expected one
    argument' before this script's own validation ever runs."""
    out = []
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok in _COORD_FLAGS and i + 1 < len(argv):
            out.append(f"{tok}={argv[i + 1]}")
            i += 2
        else:
            out.append(tok)
            i += 1
    return out


def _parse_points(s: str, n: int, flag: str) -> list:
    """Parse a flat comma-separated string of `2*n` floats into a list of `n` (x, y) tuples."""
    parts = [p.strip() for p in s.split(",")]
    if len(parts) != 2 * n:
        raise SystemExit(f"{flag} must be {2 * n} comma-separated numbers ({n} x,y pairs), "
                         f"got {len(parts)}: {s!r}")
    try:
        vals = [float(v) for v in parts]
    except ValueError:
        raise SystemExit(f"{flag} must be {2 * n} comma-separated numbers, got {s!r}")
    return [(vals[i], vals[i + 1]) for i in range(0, 2 * n, 2)]


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Automated ophys (2p) <-> Xenium coregistration for one subject.")
    ap.add_argument("--subject_id", required=True, help="Subject id, e.g. 816462")
    ap.add_argument("--input_dir", default="/root/capsule/data",
                    help="Mounted data root holding the subject's aligned-Xenium + "
                         "z-stack assets.")
    ap.add_argument("--output_dir", default="/root/capsule/results",
                    help="Output asset dir for the coregistration + QC artifacts.")
    ap.add_argument("--pose_mode", default="auto",
                    choices=["auto", "center-rotation"],
                    help="Pose-seeding mode for the anchor section (default auto = fully "
                         "automatic blind search). 'center-rotation' takes a human-provided "
                         "rough center/rotation (see --center_um/--rotation_deg/--pose_json) "
                         "and only auto-runs a depth sweep from there -- use this when 'auto' "
                         "fails to find the true pose.")
    ap.add_argument("--anchor_sec", default="",
                    help="Xenium section number to run the (expensive) initial pose search "
                         "on. Blank (default) = auto-select via "
                         "xenium_autocoreg.anchor.select_anchor_section (first cell-rich, "
                         "past-density-plateau section).")
    ap.add_argument("--center_um", default="",
                    help="REQUIRED for --pose_mode center-rotation (unless --pose_json is "
                         "given): rough anchor-section center 'X,Y' in the Xenium-aligned "
                         "frame (um), eyeballed from a confocal/vasculature image.")
    ap.add_argument("--rotation_deg", default="",
                    help="REQUIRED for --pose_mode center-rotation (unless --pose_json is "
                         "given): rough rotation (deg) of the z-stack relative to the Xenium "
                         "section.")
    ap.add_argument("--zstack_scale_to_xenium", default="",
                    help="OPTIONAL: override of the subject's own "
                         "SubjectConfig.zstack_scale_to_Xenium -- the z-stack-to-Xenium linear "
                         "PHYSICAL scale factor used to seed the initial pose search. Applies "
                         "to EVERY --pose_mode alike (not just center-rotation). NOT the same "
                         "as the fitted affine scale reported in 2p2xenium's own logs/QC (that's "
                         "a measured output, not an input). Blank (default) = use the subject's "
                         "own configured value.")
    ap.add_argument("--pose_json", default="",
                    help="OPTIONAL alternative to --center_um/--rotation_deg/"
                         "--zstack_scale_to_xenium: path to a JSON file {\"center_um\": [x,y], "
                         "\"rotation_deg\": r, \"zstack_scale_to_xenium\": s} "
                         "(\"zstack_scale_to_xenium\" optional). Keys present in the file "
                         "override the matching inline flag -- same shape "
                         "xenium_autocoreg.cli's own --pose-json takes.")
    ap.add_argument("--num_cpus", default="",
                    help="OPTIONAL: worker-process count for every parallelized stage (the auto "
                         "pose-grid search, per-section fine registration, cell-centroid "
                         "extraction, 3D point mapping). Blank/0/a value exceeding this "
                         "machine's CPU count = auto (every available core); 1 = serial, no "
                         "multiprocessing at all. See xenium_autocoreg.resources.resolve_num_cpus.")
    ap.add_argument("--zstack_registered_tif", default="",
                    help="OPTIONAL: pin the exact z-stack registered-intensity .tif to use, "
                         "bypassing subject_resolver's automatic discovery. Must be given "
                         "together with --zstack_segmented_tif (same physical acquisition, same "
                         "shape). Resolved relative to --input_dir if not absolute.")
    ap.add_argument("--zstack_segmented_tif", default="",
                    help="OPTIONAL: pin the exact z-stack segmentation-label .tif to use, "
                         "bypassing subject_resolver's automatic discovery. Must be given "
                         "together with --zstack_registered_tif. Resolved relative to "
                         "--input_dir if not absolute.")
    args = ap.parse_args(_fix_negative_coord_tokens(sys.argv[1:]))

    sid = str(args.subject_id).strip()
    if not sid:
        raise SystemExit("subject_id is required (e.g. --subject_id 816462)")
    pose_mode = args.pose_mode

    anchor_sec = None
    if args.anchor_sec:
        try:
            anchor_sec = int(args.anchor_sec)
        except ValueError:
            raise SystemExit(f"--anchor_sec must be an integer, got {args.anchor_sec!r}")
    anchor_sec_given = anchor_sec is not None

    num_cpus = None
    if args.num_cpus:
        try:
            num_cpus = int(args.num_cpus)
        except ValueError:
            raise SystemExit(f"--num_cpus must be an integer, got {args.num_cpus!r}")

    center_um, rotation_deg, zstack_scale_to_xenium = None, None, None
    if args.center_um:
        center_um = _parse_points(args.center_um, 1, "--center_um")[0]
    if args.rotation_deg:
        try:
            rotation_deg = float(args.rotation_deg)
        except ValueError:
            raise SystemExit(f"--rotation_deg must be a number, got {args.rotation_deg!r}")
    if args.zstack_scale_to_xenium:
        try:
            zstack_scale_to_xenium = float(args.zstack_scale_to_xenium)
        except ValueError:
            raise SystemExit(f"--zstack_scale_to_xenium must be a number, got "
                             f"{args.zstack_scale_to_xenium!r}")
    if args.pose_json:
        pose_json_path = Path(args.pose_json)
        if not pose_json_path.is_absolute():
            pose_json_path = Path(args.input_dir) / pose_json_path
        if not pose_json_path.exists():
            raise SystemExit(f"--pose_json not found: {args.pose_json!r} (resolved to {pose_json_path})")
        data = json.loads(pose_json_path.read_text())
        if "center_um" in data:
            center_um = tuple(data["center_um"])
        if "rotation_deg" in data:
            rotation_deg = data["rotation_deg"]
        if "zstack_scale_to_xenium" in data:
            zstack_scale_to_xenium = data["zstack_scale_to_xenium"]

    if pose_mode == "center-rotation" and (center_um is None or rotation_deg is None):
        raise SystemExit(
            "--pose_mode center-rotation requires --center_um AND --rotation_deg (or a "
            "--pose_json with both keys). Got "
            f"center_um={center_um!r} rotation_deg={rotation_deg!r}.")

    if bool(args.zstack_registered_tif) != bool(args.zstack_segmented_tif):
        raise SystemExit(
            "--zstack_registered_tif and --zstack_segmented_tif must be given together "
            f"(a mismatched auto+pinned pair would silently combine two different physical "
            f"acquisitions). Got zstack_registered_tif={args.zstack_registered_tif!r} "
            f"zstack_segmented_tif={args.zstack_segmented_tif!r}.")

    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)

    from xenium_autocoreg.anchor import select_anchor_section
    from xenium_autocoreg.cli import run_subject

    print(f"[capsule] subject={sid}  pose_mode={pose_mode}  center_um={center_um}  "
          f"rotation_deg={rotation_deg}  zstack_scale_to_xenium={zstack_scale_to_xenium}  "
          f"num_cpus={num_cpus}", flush=True)

    def _resolve_pinned(p):
        p = Path(p)
        return p if p.is_absolute() else Path(args.input_dir) / p

    zstack_pin_kwargs = {}
    if args.zstack_registered_tif and args.zstack_segmented_tif:
        reg_tif = _resolve_pinned(args.zstack_registered_tif)
        seg_tif = _resolve_pinned(args.zstack_segmented_tif)
        if not reg_tif.exists():
            raise SystemExit(f"--zstack_registered_tif not found: {args.zstack_registered_tif!r} "
                             f"(resolved to {reg_tif})")
        if not seg_tif.exists():
            raise SystemExit(f"--zstack_segmented_tif not found: {args.zstack_segmented_tif!r} "
                             f"(resolved to {seg_tif})")
        zstack_pin_kwargs = dict(zstack_registered_tif=reg_tif, zstack_segmented_tif=seg_tif)
        print(f"[capsule] PINNED z-stack override -- bypassing automatic "
              f"ophys-z-stacks_*/multiplane-ophys_* discovery ENTIRELY for this run:", flush=True)
        print(f"[capsule]   zstack_registered_tif (pinned) = {reg_tif}", flush=True)
        print(f"[capsule]   zstack_segmented_tif  (pinned) = {seg_tif}", flush=True)

    print(f"[capsule] resolving subject assets under {args.input_dir}", flush=True)
    cfg = subject_resolver.resolve_subject(int(sid), data_root=Path(args.input_dir), **zstack_pin_kwargs)

    print(f"[capsule]   aligned_dir            = {cfg.aligned_dir}", flush=True)
    print(f"[capsule]   zstack_registered_tif  = {cfg.zstack_registered_tif}", flush=True)
    print(f"[capsule]   zstack_segmented_tif   = {cfg.zstack_segmented_tif}", flush=True)
    print(f"[capsule]   zstack_xy_um           = {cfg.zstack_xy_um:.4f}", flush=True)
    print(f"[capsule]   reporter_zarr_root     = {cfg.reporter_zarr_root}", flush=True)
    print(f"[capsule]   sections ({len(cfg.sections)}) = {cfg.sections}", flush=True)

    if anchor_sec is None:
        print("[capsule] auto-selecting anchor section (xenium_autocoreg.anchor.select_anchor_section)", flush=True)
        anchor_sec, _ = select_anchor_section(cfg, verbose=True)
    print(f"[capsule] anchor section = {anchor_sec}", flush=True)

    print(f"[capsule] running xenium_autocoreg.cli.run_subject -> {out}", flush=True)
    summary = run_subject(
        cfg, out, pose_mode=pose_mode, anchor_sec=anchor_sec,
        center_um=center_um, rotation_deg=rotation_deg,
        zstack_scale_to_xenium=zstack_scale_to_xenium,
        num_cpus=num_cpus, verbose=True)

    # The package's own working tree for chain_refine (_chain_internal/) is regenerable
    # scratch, not a scientific output -- move it out of the results asset (per the
    # 2p2xenium README, which itself labels this folder "harmless scratch").
    chain_internal = out / "_chain_internal"
    if chain_internal.exists():
        scratch_dir = Path("/root/capsule/scratch/xenium_autocoreg_work") / sid
        scratch_dir.mkdir(parents=True, exist_ok=True)
        dest = scratch_dir / "_chain_internal"
        if dest.exists():
            shutil.rmtree(dest)
        shutil.move(str(chain_internal), str(dest))
        print(f"[capsule] moved regenerable working tree -> {dest}", flush=True)

    sections = [rec["sec"] for rec in summary]
    n_valid_matches = int(sum(int(rec.get("n_valid_matches", 0)) for rec in summary))
    manifest = {
        "subject_id": sid,
        "pose_mode": pose_mode,
        "anchor_sec": anchor_sec,
        "anchor_sec_auto_selected": not anchor_sec_given,
        "num_cpus": num_cpus,
        "center_um": list(center_um) if center_um else None,
        "rotation_deg": rotation_deg,
        "zstack_scale_to_xenium": zstack_scale_to_xenium,
        "zstack_pinned": bool(args.zstack_registered_tif and args.zstack_segmented_tif),
        "aligned_dir": str(cfg.aligned_dir),
        "zstack_registered_tif": str(cfg.zstack_registered_tif),
        "zstack_segmented_tif": str(cfg.zstack_segmented_tif),
        "zstack_xy_um": cfg.zstack_xy_um,
        "reporter_zarr_root": (str(cfg.reporter_zarr_root) if cfg.reporter_zarr_root else None),
        "sections": sections,
        "n_sections": len(sections),
        "n_valid_matches_total": n_valid_matches,
        "propagation_summary": summary,
    }
    (out / "coreg_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"[capsule] wrote coreg_manifest.json ({len(sections)} sections, "
          f"{n_valid_matches} total valid matches)", flush=True)

    print(f"[capsule] subject {sid} done. Outputs in {out}:", flush=True)
    for f in sorted(out.rglob("*")):
        if f.is_file():
            print("    ", f.relative_to(out), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
