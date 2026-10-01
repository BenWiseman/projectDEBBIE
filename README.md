# DEBBIE

**D**ual **E**nsemble of **B**oosted trees and **B**idirectional recurrent
networks for **I**dentifying **E**xtragalactic sources.

A catalogue of **161,721 low-surface-brightness galaxy candidates** with no prior
identification as extragalactic objects, selected from DESI Legacy Imaging
Surveys DR10 south over the complete |b| ≥ 30° footprint (≈12,900 deg²), plus a
separately flagged low-latitude extension of 78,713.

Run from public data on one desktop.

---

## The result in one line

Describing a source by its light profile in concentric annuli, and classifying
those features with gradient-boosted trees, scores **AUC 0.9810 [0.9765, 0.9852]**
on the DeepShadows benchmark, against **0.9098** for six conventional catalogue
features through the identical classifier. The gain is in the representation.
No GPU required.

## What is in the catalogue, and how clean it is

Purity is **measured**, not asserted: 999 objects from a 1,999-object sample
stratified over 28 cells of the (r_eff, μ_eff) plane, judged by eye, weighted so
the catalogue total reconstructs.

| subset | entries | purity (95% CI) |
|---|---|---|
| **μ_eff < 25.5** | **63,272** | **95.6% [92.8, 97.4]** |
| μ_eff ≥ 25.5 | 98,449 | 50.9% [45.9, 55.8] |
| all | 161,721 | 70.3% [66.8, 73.5] |

The bright subset is defined by a cut on a released column, not a model score, so
the published selection function applies to it unchanged.

191 of the 999 judged were undecided. Counting them all as contaminants gives
56.1%, all as galaxies 76.3%; the true value is inside that range.

## What the claim is, precisely

These objects are *previously unidentified as galaxies*. **99.3% carry a prior
automated designation**, most often WISE, and every one is recorded per object.
The claim is not that no instrument recorded the light. It is that no publication
or catalogue ever classified it as a galaxy. A check that *failed* is never
counted as novelty.

All redshifts are photometric. Objects are reported as **candidates**.

---

## Layout

| path | what |
|---|---|
| `catalogue/` | the catalogue, checksums, purity measurements, visual verdicts |
| `paper/` | the manuscript (Quarto), figures, bibliography |
| `code/` | the pipeline, in the order it runs |
| `showcase/` | self-contained HTML viewers for the catalogue and the screening |
| `docs/` | column definitions and methods notes |

## Verifying the catalogue

```bash
cd catalogue && sha256sum -c SHA256
```

Every number in the manuscript is filled from `catalogue/PAPER_NUMBERS.json` by
`paper/fill_numbers.py`, which **fails loudly on any unfilled placeholder**. The
paper previously quoted four different catalogue sizes at once; this makes that
impossible.

Reproduce the designation counts:

```bash
python3 - <<'PY'
import csv, io, re
def scan(p):
    return [int(re.match(r'[A-Z]+\s+(\d+)$', r['name']).group(1))
            for r in csv.DictReader(io.open(p, encoding='utf-8'))]
a, b = scan('catalogue/DEBBIE.csv'), scan('catalogue/DEBBIE_lowb_extension.csv')
print(len(a), len(b), len(a)+len(b), max(a+b), len(set(a) & set(b)))
PY
```

Expected: `161721 78713 240434 268380 0`.

## Reproducing the selection

The raw sweeps are 2.44 TB and are not distributed. `code/distil.py` produces the
3.2 GB extract (33.5M resolved sources, 31 columns) that every selection runs
against; any alternative cut re-runs from it in seconds.

## Known limitations

Photometry is the survey's and is not independently remeasured. No redshift is
spectroscopic. The selection collapses above r_eff ≈ 15″, excludes |b| < 30°
entirely, and is spatially non-uniform where prior surveys are deepest. Galactic
extinction is **not** removed before the surface-brightness cut, so a small
fraction of entries are intrinsically brighter than the stated window.

## Data credit

Any publication from this must carry the Legacy Surveys acknowledgement
verbatim. It is long and it is required.
