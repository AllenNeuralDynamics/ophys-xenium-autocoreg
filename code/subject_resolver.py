"""AIND/CodeOcean-specific subject-asset resolver for the ophys<->Xenium autocoreg capsule.

`xenium_autocoreg` (2p2xenium) is data-layout-agnostic -- it only knows about
`xenium_autocoreg.config.SubjectConfig`. This module is the CodeOcean-specific glue: it resolves
one subject's mounted data assets (attached under a capsule's /root/capsule/data, or wherever
`data_root` points) into a `SubjectConfig`, using this lab's asset-naming convention. If your data
isn't laid out this way, build a `SubjectConfig` yourself (or write your own resolver) instead of
using this module -- see the 2p2xenium README's "Configuration" section.
"""
import glob
from pathlib import Path
from typing import Optional

from xenium_autocoreg.config import SubjectConfig, zstack_xy_um_from_roi_metadata


def _latest_glob(pattern: str) -> Optional[str]:
    matches = sorted(glob.glob(pattern))
    return matches[-1] if matches else None


def _candidate_stacks(root, fov_tag=None):
    """Subdirectories of `root`, optionally filtered to ones whose name contains `fov_tag`."""
    stacks = [p for p in glob.glob(f"{root}/*") if Path(p).is_dir()]
    return sorted(s for s in stacks if fov_tag is None or fov_tag in s)


def _n_rois(seg_tif_path):
    import tifffile as tiff
    import numpy as np
    return int(np.unique(tiff.imread(seg_tif_path)).size) - 1   # exclude background (0)


def _find_zstack_pair_multiplane(data_root: Path, subject_id: int, fov_tag="700x700"):
    """Resolver for a `multiplane-ophys_{subject}_*_cortical-zstack-{segmentation,registration}_*`
    asset layout, reading pixel size from the acquisition's own metadata rather than assuming one."""
    seg_dirs = glob.glob(str(data_root / f"multiplane-ophys_{subject_id}_*_cortical-zstack-segmentation_*"))
    seg_dirs_tagged = [d for d in seg_dirs if fov_tag in d]
    seg_dirs = seg_dirs_tagged if seg_dirs_tagged else seg_dirs
    seg_tif = None
    for d in sorted(seg_dirs):
        t = _latest_glob(f"{d}/channel_0_ref_0/segmentation_masks.tif")
        if t:
            seg_tif = t
    if seg_tif is None:
        raise FileNotFoundError(f"no multiplane-ophys_{subject_id}_*_cortical-zstack-segmentation_* "
                                f"segmentation_masks.tif found under {data_root}")

    reg_dirs = glob.glob(str(data_root / f"multiplane-ophys_{subject_id}_*_cortical-zstack-registration_*"))
    reg_dirs_tagged = [d for d in reg_dirs if fov_tag in d]
    reg_dirs = reg_dirs_tagged if reg_dirs_tagged else reg_dirs
    reg_tif, roi_meta = None, None
    for d in sorted(reg_dirs):
        t = _latest_glob(f"{d}/cortical_zstack_0/channel_0_ref_0/*_2xREG.tif")
        if t:
            reg_tif = t
            roi_meta = _latest_glob(f"{d}/cortical_zstack_0/roi_groups_metadata.json")
    if reg_tif is None:
        raise FileNotFoundError(f"no multiplane-ophys_{subject_id}_*_cortical-zstack-registration_* "
                                f"*_2xREG.tif found under {data_root}")
    if roi_meta is None:
        raise FileNotFoundError(f"no roi_groups_metadata.json alongside {reg_tif} -- cannot "
                                f"determine this acquisition's um/px without it")
    return Path(reg_tif), Path(seg_tif), zstack_xy_um_from_roi_metadata(roi_meta)


def _find_zstack_pair(data_root: Path, subject_id: int, fov_tag="700x700"):
    """Picks the acquisition (matching `fov_tag` when more than one is mounted) with the most
    segmented ROIs, preferring a separate registered-intensity asset and falling back to a
    co-located raw-data file if none is mounted. Falls back to `_find_zstack_pair_multiplane` for
    a different asset-naming convention if the primary one isn't found."""
    seg_dirs = glob.glob(str(data_root / f"ophys-z-stacks_{subject_id}_segmented*"))
    if not seg_dirs:
        return _find_zstack_pair_multiplane(data_root, subject_id, fov_tag)

    seg_stacks = []
    for d in seg_dirs:
        seg_stacks += _candidate_stacks(d, fov_tag)
    if not seg_stacks:
        raise FileNotFoundError(f"no {fov_tag} segmented stack for subject {subject_id}")

    seg_tifs = {}
    for stack in seg_stacks:
        t = _latest_glob(f"{stack}/channel_0_ref_0/segmentation_masks.tif")
        if t:
            seg_tifs[stack] = t
    if not seg_tifs:
        raise FileNotFoundError(f"no segmentation_masks.tif (channel_0_ref_0) under any {fov_tag} "
                                f"stack for subject {subject_id}")

    best_stack = max(seg_tifs, key=lambda s: _n_rois(seg_tifs[s]))
    seg_tif = seg_tifs[best_stack]

    stem_name = Path(best_stack).name.split("_segmented")[0]
    reg_dir = None
    reg_dirs = glob.glob(str(data_root / f"ophys-z-stacks_{subject_id}_registered*"))
    reg_tif = None
    for d in reg_dirs:
        found_dir = _latest_glob(f"{d}/{stem_name}_registered_*")
        t = _latest_glob(f"{found_dir}/channel_0_ref_0/*_2xREG.tif") if found_dir else None
        if t:
            reg_tif, reg_dir = t, found_dir
            break
    if reg_tif is None:
        reg_tif = _latest_glob(f"{best_stack}/channel_0_ref_0/zstack_data.tif")
    if reg_tif is None:
        raise FileNotFoundError(f"no raw intensity source (registered *_2xREG.tif or co-located "
                                f"zstack_data.tif) for {best_stack}")

    # roi_groups_metadata.json's location varies by acquisition pipeline version -- check every
    # plausible spot (registered-dir root, registered-dir/channel_0_ref_0, segmented-dir/channel_0_ref_0)
    # before giving up.
    roi_meta = None
    for candidate in (
        f"{reg_dir}/roi_groups_metadata.json" if reg_dir else None,
        f"{reg_dir}/channel_0_ref_0/roi_groups_metadata.json" if reg_dir else None,
        f"{best_stack}/channel_0_ref_0/roi_groups_metadata.json",
        f"{best_stack}/roi_groups_metadata.json",
    ):
        if candidate and _latest_glob(candidate):
            roi_meta = _latest_glob(candidate)
            break
    if roi_meta is not None:
        zstack_xy_um = zstack_xy_um_from_roi_metadata(roi_meta)
    else:
        zstack_xy_um = None  # caller decides whether/how to fall back -- see resolve_subject
    return Path(reg_tif), Path(seg_tif), zstack_xy_um


def _zstack_xy_um_near(reg_tif: Path) -> Optional[float]:
    """Look for a `roi_groups_metadata.json` near an arbitrary (possibly pinned) registered-tif
    path -- checks its own directory and one level up, the same two shapes `_find_zstack_pair`
    checks for an auto-discovered stack, but relative to the file itself rather than a known
    `{stack}` folder name. Returns None (not an error) if none is found -- the caller decides
    whether/how to fall back."""
    for candidate_dir in (reg_tif.parent, reg_tif.parent.parent):
        meta = _latest_glob(str(candidate_dir / "roi_groups_metadata.json"))
        if meta:
            return zstack_xy_um_from_roi_metadata(meta)
    return None


def resolve_subject(subject_id: int, data_root: Path, fov_tag: str = "700x700",
                    fallback_fov_um: Optional[float] = 700.0,
                    fallback_native_px: Optional[int] = 512,
                    zstack_registered_tif: Optional[Path] = None,
                    zstack_segmented_tif: Optional[Path] = None,
                    zstack_xy_um: Optional[float] = None) -> SubjectConfig:
    """Resolve one subject's `SubjectConfig` from this lab's CodeOcean mounted-asset naming
    convention under `data_root` (typically /root/capsule/data).

    `fallback_fov_um`/`fallback_native_px`: used ONLY when no `roi_groups_metadata.json` can be
    found for this acquisition (some older assets don't have one mounted) -- an explicit, visible,
    overridable nominal value rather than a silent assumption. Pass `fallback_fov_um=None` to
    require real metadata and raise instead.

    `zstack_registered_tif`/`zstack_segmented_tif`: OPTIONAL explicit pin. When BOTH are given,
    this COMPLETELY BYPASSES the automatic `ophys-z-stacks_*`/`multiplane-ophys_*` discovery
    (`_find_zstack_pair`/`_find_zstack_pair_multiplane`) for this subject -- no z-stack asset needs
    to be auto-discoverable (or even attached) at all; only the pinned files are read. Use this
    for a registration/segmentation pair the automatic discovery can't reach (e.g. a
    differently-structured derived asset). `zstack_xy_um`: optional explicit calibration override
    for the pinned pair; when omitted, it's auto-derived from a `roi_groups_metadata.json` found
    near the pinned registered tif (see `_zstack_xy_um_near`), falling back to
    `fallback_fov_um`/`fallback_native_px` if none is found there either."""
    if (zstack_registered_tif is None) != (zstack_segmented_tif is None):
        raise ValueError(
            "zstack_registered_tif and zstack_segmented_tif must be given together (a "
            "mismatched pinned+auto-discovered pair would silently combine two different "
            f"physical acquisitions). Got zstack_registered_tif={zstack_registered_tif!r} "
            f"zstack_segmented_tif={zstack_segmented_tif!r}.")

    data_root = Path(data_root)
    aligned = _latest_glob(str(data_root / f"Xenium-ophys-coregistered_{subject_id}_*"))
    if aligned is None:
        raise FileNotFoundError(f"no Xenium-ophys-coregistered_{subject_id}_* aligned-frame directory "
                                f"under {data_root}")

    if zstack_registered_tif is not None and zstack_segmented_tif is not None:
        zstack_reg = Path(zstack_registered_tif)
        zstack_seg = Path(zstack_segmented_tif)
        if not zstack_reg.exists():
            raise FileNotFoundError(f"pinned zstack_registered_tif not found: {zstack_reg}")
        if not zstack_seg.exists():
            raise FileNotFoundError(f"pinned zstack_segmented_tif not found: {zstack_seg}")
        if zstack_xy_um is None:
            zstack_xy_um = _zstack_xy_um_near(zstack_reg)
            if zstack_xy_um is None:
                if fallback_fov_um is None:
                    raise FileNotFoundError(
                        f"no roi_groups_metadata.json found near pinned zstack_registered_tif="
                        f"{zstack_reg}, and no fallback_fov_um given -- pass zstack_xy_um "
                        f"explicitly")
                print(f"[{subject_id}] WARNING: no roi_groups_metadata.json found near pinned "
                      f"zstack_registered_tif -- assuming a nominal {fallback_fov_um}um FOV over "
                      f"{fallback_native_px}px (pass zstack_xy_um= to resolve_subject to override).",
                      flush=True)
                zstack_xy_um = fallback_fov_um / fallback_native_px
    else:
        zstack_reg, zstack_seg, zstack_xy_um_found = _find_zstack_pair(data_root, subject_id, fov_tag)
        if zstack_xy_um is None:
            zstack_xy_um = zstack_xy_um_found
        if zstack_xy_um is None:
            if fallback_fov_um is None:
                raise FileNotFoundError(f"no roi_groups_metadata.json found for subject {subject_id}'s "
                                        f"z-stack, and no fallback_fov_um given -- cannot determine um/px")
            print(f"[{subject_id}] WARNING: no roi_groups_metadata.json found -- assuming a nominal "
                  f"{fallback_fov_um}um FOV over {fallback_native_px}px (pass fallback_fov_um= to "
                  f"resolve_subject to change this).", flush=True)
            zstack_xy_um = fallback_fov_um / fallback_native_px

    reporter_root = _latest_glob(str(data_root / f"Xenium_{subject_id}_*_processed"))

    return SubjectConfig(
        subject_id=subject_id,
        aligned_dir=Path(aligned),
        zstack_registered_tif=zstack_reg,
        zstack_segmented_tif=zstack_seg,
        zstack_xy_um=zstack_xy_um,
        reporter_zarr_root=Path(reporter_root) if reporter_root else None,
    )


if __name__ == "__main__":
    import sys
    cfg = resolve_subject(int(sys.argv[1]), data_root=Path(sys.argv[2] if len(sys.argv) > 2 else "/root/capsule/data"))
    print(cfg)
    print(f"{len(cfg.sections)} sections: {cfg.sections}")
