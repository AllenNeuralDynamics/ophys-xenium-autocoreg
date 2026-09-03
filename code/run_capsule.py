"""Reproducible-run entry for the automated ophys (2p) <-> Xenium coregistration capsule.

Wraps the `xenium_autocoreg` package from
https://github.com/jkim0731/2p2xenium (installed in environment/postInstall) for one
subject. The whole pipeline is GT-free: anchor-section selection, automatic pose
seeding (soma-print point-cloud descriptor), 3D tilt fitting, chain-refine propagation
across sections, fine mask-TPS registration + probability-filtered cell matching, and
3D point mapping of every Xenium cell into z-stack coordinates.

Pipeline stages (see the 2p2xenium README for full detail; all run via
``xenium_autocoreg.cli.run_subject``):
    1. anchor-section selection      2. pose seeding (auto / center-rotation / corners)
    3. tilt fitting                  4. chain-refine propagation (all sections)
    5. fine registration + cell matching (per section)
    6. 3D point mapping of every Xenium cell centroid

Inputs (attach as data assets, mounted under ``--input_dir`` = /root/capsule/data; resolved
by ``xenium_autocoreg.config.resolve_subject`` -- glob, no IDs hard-coded):
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
"""
import argparse
import json
import os
import shutil
import sys
from pathlib import Path


def _parse_center_um(s: str) -> tuple:
    parts = [p.strip() for p in s.split(",")]
    if len(parts) != 2:
        raise SystemExit(f"--center_um must be 'X,Y' (comma-separated), got {s!r}")
    try:
        return (float(parts[0]), float(parts[1]))
    except ValueError:
        raise SystemExit(f"--center_um must be two numbers 'X,Y', got {s!r}")


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
                    choices=["auto", "center-rotation", "corners"],
                    help="Pose-seeding mode for the anchor section (default auto = fully "
                         "automatic blind search). 'center-rotation' takes a human-provided "
                         "rough center/rotation (see --center_um/--rotation_deg/--pose_json) "
                         "and only auto-runs a depth sweep from there -- use this when 'auto' "
                         "fails to find the true pose. 'corners' is NOT YET IMPLEMENTED "
                         "upstream (2p2xenium xenium_autocoreg.pose_seed.seed_from_corners).")
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
                         "override the matching inline flag.")
    args = ap.parse_args()

    sid = str(args.subject_id).strip()
    if not sid:
        raise SystemExit("subject_id is required (e.g. --subject_id 816462)")
    pose_mode = args.pose_mode

    if pose_mode == "corners":
        raise SystemExit(
            "--pose_mode corners is not yet implemented upstream -- see "
            "xenium_autocoreg.pose_seed.seed_from_corners in 2p2xenium for exactly why, "
            "and what real click data is needed before it can be. Use 'auto' or "
            "'center-rotation' instead.")

    anchor_sec = None
    if args.anchor_sec:
        try:
            anchor_sec = int(args.anchor_sec)
        except ValueError:
            raise SystemExit(f"--anchor_sec must be an integer, got {args.anchor_sec!r}")
    anchor_sec_given = anchor_sec is not None

    center_um, rotation_deg, scale = None, None, None
    if args.center_um:
        center_um = _parse_center_um(args.center_um)
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

    # Configure the package via env BEFORE importing it (config.py reads
    # XENIUM_AUTOCOREG_DATA_ROOT at import time to build its DATA_ROOT constant).
    os.environ["XENIUM_AUTOCOREG_DATA_ROOT"] = str(args.input_dir)

    from xenium_autocoreg.config import resolve_subject
    from xenium_autocoreg.anchor import select_anchor_section
    from xenium_autocoreg.cli import run_subject

    print(f"[capsule] subject={sid}  pose_mode={pose_mode}  "
          f"center_um={center_um}  rotation_deg={rotation_deg}  scale={scale}", flush=True)

    print(f"[capsule] resolving subject assets under {args.input_dir}", flush=True)
    cfg = resolve_subject(int(sid))
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
        int(sid), out, pose_mode=pose_mode, anchor_sec=anchor_sec,
        center_um=center_um, rotation_deg=rotation_deg, scale=scale, verbose=True)

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
