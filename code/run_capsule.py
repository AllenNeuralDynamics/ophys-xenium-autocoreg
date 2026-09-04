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
    center-rotation  -- ``--center_um`` + ``--rotation_deg`` required (``--scale`` optional).

(A third, corner-based mode is a documented TODO in ``xenium_autocoreg.pose_seed.seed_from_corners``
-- not yet implemented upstream, and not exposed here.)
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
    ap.add_argument("--scale", default="",
                    help="OPTIONAL for --pose_mode center-rotation: z-stack-to-Xenium scale "
                         "factor. Blank (default) = the subject's own "
                         "SubjectConfig.zstack_scale_to_Xenium.")
    ap.add_argument("--pose_json", default="",
                    help="OPTIONAL alternative to --center_um/--rotation_deg/--scale: path "
                         "to a JSON file {\"center_um\": [x,y], \"rotation_deg\": r, "
                         "\"scale\": s} (\"scale\" optional). Keys present in the file "
                         "override the matching inline flag -- same shape "
                         "xenium_autocoreg.cli's own --pose-json takes.")
    ap.add_argument("--num_cpus", default="",
                    help="OPTIONAL: worker-process count for every parallelized stage (the auto "
                         "pose-grid search, per-section fine registration, cell-centroid "
                         "extraction, 3D point mapping). Blank/0/a value exceeding this "
                         "machine's CPU count = auto (every available core); 1 = serial, no "
                         "multiprocessing at all. See xenium_autocoreg.resources.resolve_num_cpus.")
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

    center_um, rotation_deg, scale = None, None, None
    if args.center_um:
        center_um = _parse_points(args.center_um, 1, "--center_um")[0]
    if args.rotation_deg:
        try:
            rotation_deg = float(args.rotation_deg)
        except ValueError:
            raise SystemExit(f"--rotation_deg must be a number, got {args.rotation_deg!r}")
    if args.scale:
        try:
            scale = float(args.scale)
        except ValueError:
            raise SystemExit(f"--scale must be a number, got {args.scale!r}")
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
        if "scale" in data:
            scale = data["scale"]

    if pose_mode == "center-rotation" and (center_um is None or rotation_deg is None):
        raise SystemExit(
            "--pose_mode center-rotation requires --center_um AND --rotation_deg (or a "
            "--pose_json with both keys). Got "
            f"center_um={center_um!r} rotation_deg={rotation_deg!r}.")

    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)

    from xenium_autocoreg.anchor import select_anchor_section
    from xenium_autocoreg.cli import run_subject

    print(f"[capsule] subject={sid}  pose_mode={pose_mode}  center_um={center_um}  "
          f"rotation_deg={rotation_deg}  scale={scale}  num_cpus={num_cpus}", flush=True)

    print(f"[capsule] resolving subject assets under {args.input_dir}", flush=True)
    cfg = subject_resolver.resolve_subject(int(sid), data_root=Path(args.input_dir))
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
        center_um=center_um, rotation_deg=rotation_deg, scale=scale,
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
        "scale": scale,
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
