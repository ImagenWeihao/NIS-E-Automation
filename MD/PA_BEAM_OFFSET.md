# 850 nm PA beam offset — measured 2026-08-07, sph_A02_01 (well A02)

> **Status: MEASURED, NOT CALIBRATED. Do not hard-code these numbers into the PA job
> or the GUI yet.** One spheroid, one run, and the Z figure is confounded (see
> [Confounds](#confounds)). This document exists so the correction is implemented from
> evidence rather than from memory, and so the follow-up run is designed to remove the
> confound.

The 850 nm activation optical config (`850 nm power loop full reso2`) does **not**
photoactivate where the job is told to. Activation landed **below and to one side** of
the commanded PA volume.

## Commanded vs observed

Commanded values are read from the PA stack's own metadata, not from the GUI fields.

| | Commanded | Observed | Offset |
|---|---|---|---|
| Z (centre) | 7849.99 µm | **7813.9 µm** | **−36.1 µm** |
| Z (bottom plane) | 7789.98 µm | — | activation extends to 7780, i.e. ~10 µm *below* the lowest commanded plane |
| X | 39752.56 µm (square centre) | — | **−0.3 µm** (on target) |
| Y | 29095.84 µm (square centre) | — | **−47.7 µm** |

PA volume as commanded: 79.5 µm square × 25 planes, 7789.98–7910.00 µm, 5 µm step,
Bottom-to-Top, 70 loops/plane.

**The Y offset (47.7 µm) exceeds the square's half-width (39.8 µm)** — so the activated
population sits *beyond* the square edge on one side, not distributed around it. This is
a directional beam offset, not "the T-cells happen to be peripheral".

## Evidence A — activation, 940 nm PAsfGFP (post − pre)

Blobs = `post − pre > 80`, 3×3 binary opening, connected-component labelled. Signed
distance is to the PA square edge; negative = inside.

| z (µm) | blobs | px | peak | edge dist µm (mean / min / max) |
|---|---|---|---|---|
| 7780 | 1 | 24 | 222 | −11.5 / −13.1 / −9.9 |
| 7785 | 2 | 57 | 1310 | −10.4 / −13.7 / −3.7 |
| 7790 | 1 | 45 | 150 | −4.0 / −6.2 / −1.9 |
| 7795 | 2 | 18 | 164 | +13.4 / +8.1 / +18.6 |
| 7800 | 2 | 24 | 533 | +18.2 / +17.4 / +19.3 |
| 7815 | 1 | 22 | 255 | +27.7 / +26.7 / +28.6 |
| 7820 | 2 | **174** | 739 | +34.5 / +31.1 / +39.8 |
| 7825 | 5 | **182** | 527 | +37.5 / +31.7 / +51.6 |
| 7835 | 2 | 35 | 167 | +21.8 / +19.9 / +23.6 |
| 7840 | 1 | 14 | 131 | +27.7 / +26.7 / +28.6 |

Band 7780–7840 µm; intensity-weighted centroid **7813.9 µm**; weighted mean edge
distance **+23.9 µm (outside)**. Every other plane of 49 has post ≤ pre.

## Evidence B — beam footprint, 1050 nm bleach imprint

Independent of any activation signal: mean post/pre inside the square vs a matched ring
just outside it, over structure pixels only (1050 nm > p90).

| z (µm) | square | ring | imprint (sq/ring) |
|---|---|---|---|
| 7780 | 0.436 | 0.443 | 0.98 |
| 7785 | 0.224 | 0.768 | **0.29** |
| **7790** | 0.095 | 0.419 | **0.23** ← sharpest |
| 7795 | 0.112 | 0.294 | 0.38 |
| 7800 | 0.126 | 0.357 | 0.35 |
| 7810 | 0.230 | 0.460 | 0.50 |
| 7850 | 0.478 | 0.550 | 0.87 (commanded centre — imprint nearly gone) |

The square is bleached up to **4.4× harder** than immediately outside it, peaking at
7790 and fading by 7850. This is the visible "rectangle indentation" in the depth
montages.

## Efficiency

Whole-spheroid 940 nm post/pre over the 1050 nm structure mask (NOT in-square vs
out-of-square — see [Method notes](#method-notes)):

| | planes | 940 post/pre |
|---|---|---|
| activation band 7780–7840 | 13 | 0.855 |
| control ≥ 7880 | 19 | 0.949 |
| **PA efficiency (band / control)** | | **0.90×** |

Below 1.0: across the whole spheroid, bleaching outweighs conversion. Conversion is real
but confined to a few blobs, while 25 planes × 70 loops bleached the square to 0.10–0.23
of baseline.

## Confounds

1. **Z offset vs acquisition order.** The stack ran Bottom-to-Top, so the first plane
   (7790) always meets the freshest fluorophore and bleaches hardest. The footprint
   peaking at 7790 is therefore *partly* an ordering artifact. It cannot currently be
   separated from a genuine focal offset.
2. **25 planes cannot distinguish** "the focus is offset" from "only the first plane or
   two converted, the rest bleached what was already made".
3. **n = 1 spheroid, 1 depth, 1 well.** No evidence yet that the offset is constant
   across the plate or with depth.
4. The Y offset is measured from activation-blob centroids, which are biased toward
   wherever PAsfGFP-expressing cells actually are. It is strongly suggestive (47.7 µm in
   Y vs 0.3 µm in X) but should be confirmed on a target that isn't cell-distribution
   dependent — the bleach footprint is the better probe.

## Proposed calibration run (do this before the batch)

Single run, answers Z and XY at once, removes confound 1 and 2:

1. **Narrow PA: 3 planes, not 25**, at a known commanded z. One plane is enough for the
   footprint; three gives the band shape.
2. Keep 70 loops so the per-plane dose matches production.
3. **Post-PA z-stack much wider than the PA range** (±120 µm as now) so the activation
   band can be *located*, not assumed.
4. Measure both:
   - activation-blob centroid (940 post − pre) → XY and Z of the converted volume
   - bleach-square footprint (1050 post/pre, square vs ring) → XY and Z of the beam,
     independent of cell distribution
5. Repeat on a second spheroid at a different depth before generalising.

Only then fold the offset into the job/GUI.

## Method notes

- **Do not score a single "middle" plane.** The effect is ~10–60 µm from the commanded
  centre; scoring 7850 alone shows nothing.
- **Do not average over the whole PA z-range.** A ~10 µm band diluted across 120 µm
  vanishes. Localise in depth first, then measure.
- **Do not use "outside the square" as the control.** The activated cells are outside the
  square edge, so that region contains signal. Using it put signal in the denominator and
  produced selectivity 0.55 — an artifact. Use planes the beam never reached instead.
- **Connected-component label before dismissing a peak.** The 1311-count peak at 7785
  looked like a hot pixel; it is a 57-px connected blob.
- **Fixed, identical display range for pre and post.** A percentile stretch renormalises
  each panel to its own noise floor and manufactures structure in empty channels.

## Data and reproduction

| | |
|---|---|
| PA stack (25 planes) | `work/PA/20260807_155646_016__Point0000_ZStack00NN_*.nd2` |
| Pre-PA validation | `work/0807/nd2/prePA_{890nm_mBeRFP,940nm_PAsfGFP,1050nm_depth}/` |
| Post-PA validation | `work/0807/nd2/postPA_{...}/` |
| Depth montages (native 1:1) | `work/0807/sph_A02_01_depth_montage_*_native.png` |
| Duplicates quarantined | `work/0807/_recycle/` (byte-identical copies of the postPA stacks, written by a stale ND Save-to-File path) |

Acquisition conditions were verified identical pre vs post before any comparison:
same 49×512×512 geometry, 0.62148 µm/px, 318.2 µm FOV, z 7730–7970 @5 home 7850, same
stage XY, same laser powers (890/940 at 30 %, 1050 at 50 %).

Analysis scripts (scratchpad, 2026-08-07): `pa_spatial.py` (footprint, blobs, XY offset,
efficiency), `pa_depth_montage_native.py` (montages), `pa_final.py` (efficiency figure).

## Related

- `nise-850pa-oc-beam-offset` (memory)
- `nise-pa-point-from-ndtopointset` — the PA point comes from `PredefinedPoints`, not the
  ND point set
- `nise-ir-shutter-not-automatable`
