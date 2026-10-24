# DEBBIE

**D**R10 **E**nsemble of **B**oosted trees and **B**idirectional recurrent networks for
**I**dentifying **E**xtragalactic sources.

DEBBIE v1.5 holds 427,022 low surface brightness galaxy candidates from the DESI Legacy
Imaging Surveys DR10 south sweeps that no galaxy catalogue I searched, NED, SIMBAD or DESI DR1
identifies as a galaxy. Designations, DEBBIE NNNNNN, are in the IAU Dictionary of Nomenclature
(accepted 2026-10-22).

| part | Galactic latitude | entries | purity, % (95% interval) |
|---|---|---|---|
| original | \|b\| >= 30 | 161,533 | 89.1 (86.4–91.6) |
| recovered | \|b\| >= 30 | 141,405 | 86.1 (81.1–89.9) |
| extension | 12 < \|b\| < 30 | 124,084 | 84.8 (80.3–88.4) |

Purity is the weighted share of visually inspected objects with a galaxy at the catalogued
position, undecided objects left out. The recovered and extension figures come from 800 of a
2,000-object sample, judged in random order.

## Submitting

`SUBMIT.md` has the Zenodo upload, the Research Note submission with its cover letter, and
what follows publication.

## Layout

| path | what |
|---|---|
| `catalogue/` | v1.5: `DEBBIE_v1.5.csv.gz`, its `README.md` (the Zenodo description), `SHA256SUMS`, `checks.json` |
| `catalogue/purity/` | purity results and every visual verdict behind them |
| `catalogue/v1.2/` | the v1.2 release of 2026-10-01, superseded |
| `rnaas/` | the Research Note: `debbie_rnaas.pdf` and `.tex` to submit, the template and `fill_rnaas.py` that fills it |
| `paper/` | the AJ draft (Quarto), figures and bibliography; its numbers are still v1.2 (161,721 entries) |
| `tools/` | `debbie_images.py`: choose entries and download their DR10 image, Tractor model and residual, with a gallery page |
| `code/` | the v1.2 pipeline; the v1.5 scripts are not yet copied here |
| `showcase/` | v1.2 HTML viewers |
| `private/` | gitignored: CDS guide, venue notes, review correspondence, handover |

The raw `DEBBIE_v1.5.csv` and `DEBBIE_v1.5.ecsv` (116 and 147 MB) are over GitHub's 100 MB
file limit. They sit in `catalogue/` for the Zenodo upload and are gitignored.

## Getting images

`tools/debbie_images.py` picks entries by name, position, part or any column, can take a
random or stratified sample, and downloads each one's DR10 image, Tractor model and residual
from the Legacy Surveys cutout service, with `index.html` to look through them. Standard
library only.

```bash
# 60 entries spread over score and size, with a wider context image
python3 tools/debbie_images.py --sample 60 --stratify deepshadows_mean:4 shape_r:3 --context
# the 40 highest-scoring entries above 12"
python3 tools/debbie_images.py --where "shape_r>=12" --sort deepshadows_mean --desc --limit 40
# everything within 20' of a position, as FITS in g, r and z
python3 tools/debbie_images.py --cone 185.45 -8.85 20 --format fits
```

`python3 tools/debbie_images.py -h` lists every option. Output goes to `debbie_images/`
(gitignored). A stratified sample's `selection.csv` gives each entry a weight, the number of
selected entries it stands for. Tests: `python3 -m unittest tools/test_debbie_images.py`.

## Checking the files

```bash
cd catalogue
zcat DEBBIE_v1.5.csv.gz > DEBBIE_v1.5.csv
sha256sum -c --ignore-missing SHA256SUMS
```

`pandas.read_csv("catalogue/DEBBIE_v1.5.csv.gz")` reads the table as it is.

## Updating this repo

Work happens in `~/star`. `~/star/sweeps/sync_projectdebbie.py` copies the release, the
Note, the paper and the purity files here and checks every copy by SHA-256. It refuses to
overwrite a file with uncommitted edits here, so edit those files in `~/star`, not here.
README.md, SUBMIT.md and `tools/` are written here and the sync never touches them.

## Caveats

Photometry is the survey's. No redshift is spectroscopic. The cuts keep most known galaxies
with `SHAPE_R` of 3 to 15" in the surface-brightness window and almost none larger. The screen
is uneven across the sky: inside the SDSS footprint far fewer entries survive than outside it.
Extinction is not removed before the surface-brightness cut. Numbers are never reused, so the
sequence has gaps.

## Data credit

Work using this catalogue uses Legacy Surveys DR10, and should carry the Legacy Surveys
acknowledgement (https://www.legacysurvey.org/acknowledgment/). The screen also used DESI DR1,
DES DR2, SDSS-III and 2MASS; their texts are linked in `SUBMIT.md`.
