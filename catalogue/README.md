# DEBBIE v1.5

Low surface brightness galaxy candidates with no prior identification as galaxies, selected from
the DESI Legacy Imaging Surveys DR10 south sweep catalogues. 427,022 entries, designated
DEBBIE NNNNNN (IAU Dictionary of Nomenclature of Celestial Objects, accepted 2026-10-22).

Wiseman, B., in preparation: a Research Note of the AAS, then a full paper planned for AJ.

## Files

- `DEBBIE_v1.5.csv`: one row per entry, ordered by number.
- `DEBBIE_v1.5.ecsv`: the same table at full precision, with units and column descriptions.
- `SHA256SUMS`: checksums of the two tables.

## Selection

Every entry is a DR10 Tractor source with `TYPE` EXP, DEV, SER or REX; `SHAPE_R` >= 3";
mean g surface brightness within `SHAPE_R`, 22.5 - 2.5 log10(`FLUX_G`) + 2.5 log10(2 pi
`SHAPE_R`^2), between 24.0 and 26.5 mag arcsec^-2 with no extinction correction; `NOBS_G`
and `NOBS_R` >= 2; `FLUX_G`, `FLUX_R` and `FLUX_Z` > 0; and `MASKBITS` = 0. Each entry's
row in this table passes every one of these cuts.

An entry is removed if a catalogued galaxy (HyperLEDA, SDSS DR12, GLADE+, 2MASS XSC, 6dFGS,
SMUDGes) lies within the larger of 2 `SHAPE_R` and 10"; if it lies inside half a HyperLEDA
galaxy's D25 diameter; if NED or SIMBAD lists a galaxy, quasar or AGN within the larger of 5"
and `SHAPE_R`; if a DESI DR1 spectrum with ZWARN = 0 lies within 2"; or if three or more other
candidates lie within 1.5'. The `original` part was screened with fixed radii instead of the
NED/SIMBAD rule (NED galaxy within 30", any SIMBAD object within 15", any index object
including DES DR2 within 45").

## Parts

| part | Galactic latitude | entries |
|---|---|---|
| original | \|b\| >= 30 | 161,533 |
| recovered | \|b\| >= 30 | 141,405 |
| extension | 12 < \|b\| < 30 | 124,084 |

## Columns

`name`, `number`; `ra`, `dec` (Tractor, Gaia frame) and `position` (sexagesimal); `part`,
`gal_b`; `brickid`, `objid` (to join back to DR10); `type`, `shape_r`, `shape_e1`, `shape_e2`,
`sersic`; `mu_eff_g`; `flux_g/r/z/w1/w2` and `flux_ivar_g/r/z/w1` (nanomaggies);
`nobs_g/r/z`; `ebv` (not applied anywhere).

Scores from radial light profiles: `deepshadows_trees`, `deepshadows_network` and their mean
`deepshadows_mean`, trained on the DeepShadows labels. `stage2_*` come from a second pair that
fails against visual verdicts and should not be used for selection. No entry is removed on
its score; `nan` where the entry has no usable image.

`centroid_offset`: distance from the catalogued position to a flux-weighted centroid on the
DR10 images. `centroid_flag`: 0 measured, 1 beyond 2.5" or unmeasurable, 2 no usable image.

## Caveats

Purity from visual inspection of stratified samples is 85 to 89 per cent depending on the
part; the recovered and extension figures come from a sample still being judged. The cuts
keep most known galaxies with `SHAPE_R` of 3 to 15" in the surface-brightness window and
almost none larger. The screen is uneven across the sky: inside the SDSS footprint far fewer
entries survive than outside it. Numbers are assigned once and never reused, so the sequence
has gaps.
