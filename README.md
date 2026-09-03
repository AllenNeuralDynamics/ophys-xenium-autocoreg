# ophys <-> Xenium automated coregistration (capsule)

Reproducible-run CodeOcean capsule that **automatically coregisters** one subject's
in-vivo 2-photon (ophys) cortical z-stack to its Xenium spatial-transcriptomics sections,
so the same physical cell can be found in both modalities.

It wraps the `xenium_autocoreg` package from
[`2p2xenium`](https://github.com/jkim0731/2p2xenium) (installed in `environment/postInstall`).
The whole pipeline is **GT-free** (no manual ground-truth landmarks are consumed at run
time -- the initial pose is found automatically, or from a rough human-provided center/
rotation, or from clicked tissue corners).

`xenium_autocoreg` itself is **data-layout-agnostic** -- it only knows about
`xenium_autocoreg.config.SubjectConfig` (see its README's "Configuration" section). This
capsule's own [`code/subject_resolver.py`](code/subject_resolver.py) is the
CodeOcean/company-specific glue that resolves one subject's mounted assets (this lab's
`Xenium-ophys-coregistered_*`/`ophys-z-stacks_*` asset-naming convention) into a
`SubjectConfig` -- that resolver logic does **not** live in `xenium_autocoreg` (it moved
here specifically so the general package stays reusable outside this lab's CodeOcean layout).

> **Note:** `environment/postInstall` currently pins `git clone -b general-config` -- this
> capsule's parameter contract (see below) matches that branch's `xenium_autocoreg.cli`, not
> yet 2p2xenium's default branch. Drop `-b general-config` once that branch is merged upstream.

## Reproducible Run

Clicking **Reproducible Run** executes `code/run` -> `code/run_capsule.py` for one
subject (`--subject_id`, set via the app panel):

```
anchor-section selection -> pose seeding (auto / center-rotation / corners)
      |
tilt fitting  ->  chain-refine propagation (all sections)
      |
fine registration (mask-TPS + probability-filtered cell matching), per section
      |
3D point mapping of every Xenium cell -> z-stack coordinates
      |
coreg_manifest.json  (run parameters + per-section summary)
```

All stages run via `xenium_autocoreg.cli.run_subject` -- see the
[2p2xenium README](https://github.com/jkim0731/2p2xenium#the-process) for the full
per-stage algorithm detail.

### Parameters (`.codeocean/app-panel.json`)

`pose_mode` selects which of the 3 pose-seeding protocols seeds the anchor section's
initial pose. Parameter names below match `xenium_autocoreg.cli`'s own CLI/`--pose-json`
shape **exactly** (see the 2p2xenium README) -- a single `pose_json` file works unchanged
against either tool.

| param | meaning |
|---|---|
| `subject_id` | subject id, e.g. `816462` (**required**) |
| `pose_mode` | `auto` (default) fully automatic pose search / `center-rotation` human-seeded / `corners` (not yet implemented upstream, see below) |
| `anchor_sec` | Xenium section to run the initial pose search on; blank (default) = auto-select |
| `center_um` | `center-rotation` only (**required** unless `pose_json` set): rough anchor center `X,Y` (um) |
| `rotation_deg` | `center-rotation` only (**required** unless `pose_json` set): rough rotation (deg) |
| `scale` | `center-rotation`/`corners` only (optional): z-stack-to-Xenium scale factor; blank = subject default |
| `xenium_trapezoid_corners_um` | `corners` only (**required** unless `pose_json` set): 4 Xenium tissue-trapezoid corners, `X1,Y1,X2,Y2,X3,Y3,X4,Y4` (um) |
| `top_edge` | `corners` only (**required** unless `pose_json` set): which corner/edge (0-3) is the trapezoid's short/slanted top |
| `zstack_corners_um` | `corners` only (optional): the z-stack's own 4 corners, same shape as `xenium_trapezoid_corners_um`; blank = the z-stack's canonical FOV rectangle |
| `pose_json` | optional alternative to the inline fields above -- for center-rotation: `{"center_um":[x,y],"rotation_deg":r,"scale":s}`; for corners: `{"xenium_trapezoid_corners_um":[[x,y]x4],"top_edge":0,"zstack_corners_um":[[x,y]x4],"scale":s}` |

Use `center-rotation` (a rough center/rotation eyeballed from a confocal or vasculature
image) whenever `auto` fails to find the true pose -- see "Known limitations" below.
`corners` is not yet implemented upstream (`xenium_autocoreg.pose_seed.seed_from_corners`
raises `NotImplementedError`); its parameters are validated and passed through end-to-end
today (fails fast, before any expensive registration work) so no capsule changes will be
needed once it lands.

### Inputs (attach as data assets, mounted under `/root/capsule/data`)
The subject is resolved by glob via this capsule's own `subject_resolver.resolve_subject`
(no IDs hard-coded -- attach per subject, or trigger via a pipeline monitor):
- `Xenium-ophys-coregistered_{subject_id}_*` -- **required**: the aligned-frame Xenium
  directory (`Xenium_images/section_N_Neurons_aligned.tif`,
  `Xenium_segmentation_masks/section_N_Masks_aligned.tif`; sections already
  section-to-section aligned by an upstream process -- see the 2p2xenium README's
  "Required input data" for exactly what "aligned" means here).
- `ophys-z-stacks_{subject_id}_registered_*` + `ophys-z-stacks_{subject_id}_segmented_*`
  -- **required**: the registered intensity z-stack and its matching segmentation label
  volume (same shape, same physical FOV). Falls back to the
  `multiplane-ophys_{subject_id}_*cortical-zstack-{registration,segmentation}_*` naming
  if the primary one isn't attached.
- `Xenium_{subject_id}_*_processed` -- **optional**: the Xenium processed asset, used as
  the reporter-transcript (SYFP2/EGFP) population source to restrict matching to
  reporter+ cells; falls back to all segmented Xenium cells if not attached.

### Outputs (to `/root/capsule/results`, captured as a derived data asset)

Written directly by `xenium_autocoreg.cli.run_subject` in the layout documented in the
[2p2xenium README](https://github.com/jkim0731/2p2xenium#output-structure):

```
Affine matrices/        section_N_affine_matrix.npy (Xenium -> z-stack) + section_N_rotation_3d.npy
Xenium_affine_transformed/   Xenium intensity/masks, affine-warped into the z-stack canvas
warped_zstacks/          z-stack intensity + masks, TPS-warped into the shared frame
post_affine_warping/     fine tile-correlation landmark points, per section
cell_matching_probability/   per-pair IoU + Mahalanobis-distance match probabilities
cell_centroids/          Xenium + z-stack cell centroids, per section
mapped_3d_coordinates/   every Xenium cell -> 3D z-stack coordinates (non-rigid + rigid)
ophys-z-stacks/, ophys-z-stacks_segmentation_masks/   raw z-stack volume(s), copied as-is
QC/                      per-section + anchor-only QC figures (see below)
propagation_summary.json    per-section run summary (written by the package itself)
coreg_manifest.json      THIS capsule's own summary: run parameters + per-section counts
```

**Regenerable working caches**: `xenium_autocoreg`'s own `_chain_internal/` working tree
(chain-refine's internal scratch, per its own README) is moved to
`/root/capsule/scratch/xenium_autocoreg_work/<sid>/` after the run -- it is not a
scientific output.

## QC figures
- `QC/xenium_affine_zstack/section_N_xenium_affine_zstack.png` -- one per section: z-stack
  / Xenium / overlap (intensity + masks) after the final registration.
- `QC/cell_matching/section_N_cell_matching.png` -- one per section: cell contours +
  overlay of the statistically-valid matched pairs.
- `QC/initial_match/section_N_initial_match_wide.png` -- **anchor section only**: the
  initial landmark search vs. the final tilt+TPS-warped registration.

## Notes
- Part of the ophys<->Xenium coregistration workflow: this capsule produces the automated
  cell-cell matches + 3D point mapping; downstream analysis consumes
  `cell_matching_probability/` and `mapped_3d_coordinates/` directly (no separate
  interactive QC capsule is required, unlike the 2p<->3D-FISH workflow).
- This mirrors the [`capsule-2p-3DFISH-autocoreg`](https://github.com/AllenNeuralDynamics/capsule-2p-3DFISH-autocoreg)
  capsule pattern (a thin `run_capsule.py` wrapping a pip-installable `autocoreg`
  package cloned in `environment/postInstall`), applied to
  [`2p2xenium`](https://github.com/jkim0731/2p2xenium) instead of
  [`2p2fish`](https://github.com/jkim0731/2p2fish) -- with one further split: 2p2xenium
  itself carries no lab-specific asset-naming logic (unlike 2p2fish, which the 3D-FISH
  capsule wraps as-is), so this capsule owns that resolver (`code/subject_resolver.py`)
  the way `run_capsule.py` in a from-scratch capsule would.

## Known limitations
(inherited from [2p2xenium](https://github.com/jkim0731/2p2xenium#known-limitations))
- The SVD scale-clip in chain-refine bounds scale but not shear -- a weak-correlation
  section can still show shear-driven cell-shape distortion.
- `pose_mode=auto` can fail outright if the true pose sits outside the searched position
  window; switch to `pose_mode=center-rotation` when it does.
- `pose_mode=corners` is not implemented upstream (`pose_seed.seed_from_corners`) and will
  fail the run with `NotImplementedError` -- its parameters (`xenium_trapezoid_corners_um`,
  `top_edge`, `zstack_corners_um`) are already wired through end-to-end.
