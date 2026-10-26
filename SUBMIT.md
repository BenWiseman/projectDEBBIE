# Submitting DEBBIE v1.5

Order: reserve the Zenodo DOI, put it in the Note, publish Zenodo, submit the Note.

## 1. Zenodo (Ben's account)

New upload in the AAS Journals community: https://zenodo.org/uploads/new?community=aas
(the RNAAS guidelines point there for data tables).

- Files, all from `catalogue/`: `DEBBIE_v1.5.csv`, `DEBBIE_v1.5.ecsv`, `README.md`,
  `SHA256SUMS` (238 MB together); and `rnaas/debbie_rnaas.pdf`, added after the DOI is in it.
- DOI: "Do you already have a DOI?" No, then "Get a DOI now!". Send the DOI to Claude before
  publishing. It goes into `~/star/paper/rnaas/fill_rnaas.py` (ZENODO), the Note is rebuilt
  and the sync copies it here.
- Resource type: Dataset.
- Title: DEBBIE v1.5: 427,022 Galaxy Candidates with No Prior Identification in the DESI
  Legacy Imaging Surveys DR10 (the Note's title, with the version)
- Creator: Wiseman, Benjamin (ORCID 0009-0002-1023-9026). Affiliation: Linnet Labs PTY.
- Description: paste `catalogue/README.md`.
- Version: 1.5. Licence: Creative Commons Attribution 4.0 International.
- Keywords: low surface brightness galaxies; galaxy catalogues; DESI Legacy Imaging Surveys;
  DR10; DEBBIE
- Related works, once the Note has a DOI: "is supplemented by" the Note.

## 2. Research Note

Submit at https://journals.aas.org/submission/ (journal: Research Notes of the AAS).

- Files: `rnaas/debbie_rnaas.pdf` and `rnaas/debbie_rnaas.tex`. The table is in the .tex;
  there is no figure.
- Title, abstract and author come from the .tex. UAT keywords are collected at acceptance,
  not now (AAS manuscript preparation page, read 2026-10-24).
- Words: 1,424 of 1,500 by texcount with references; abstract 148 of 150. The AAS counts
  acknowledgements ("these count against your 1500 word limit", RNAAS preparation
  guidelines).

### Cover letter

Dear Editor,

I submit "DEBBIE: 427,022 Galaxy Candidates with No Prior Identification in the DESI Legacy
Imaging Surveys DR10" as a Research Note. Prof. Lintott confirmed on 20 October that a
catalogue may appear first in RNAAS and then in a main journal, and asked whether the Note
could describe the selection well enough for the catalogue to be used before the full paper.
It gives every cut and screening rule, and the catalogue, with each entry's Tractor
photometry, scores and centroid offset, is at Zenodo (doi:10.5281/zenodo.23659669). The designations are registered
in the IAU Dictionary of Nomenclature, accepted on 22 October.

The Note uses data from the Legacy Surveys, DESI DR1, DES DR2 and SDSS-III. Their required
acknowledgements come to about 1,300 words, so the Note cites each survey and links its full
text. If you need any of them in full, I will cut the body to fit.

Best wishes,
Ben Wiseman

## 3. After publication

- Dictionary revision: N = 427,022 and the Note as reference, through
  http://cdsarc.u-strasbg.fr/viz-bin/DicForm?AcA=DEBBIE&AcM=Benjamin.h.wiseman%40gmail.com
- Zenodo: add the Note's DOI under related works.

## Acknowledgement lengths (words, counted 2026-10-24)

| survey | words | source |
|---|---|---|
| DESI DR1 | 202 | https://data.desi.lbl.gov/doc/acknowledgments/ |
| Legacy Surveys, full | 652 | https://www.legacysurvey.org/acknowledgment/ |
| Legacy Surveys, short (headed for DESI collaborators) | 230 | same |
| DES | 311 | https://www.darkenergysurvey.org/the-des-project/data-access/ |
| SDSS-III | 147 | https://www.sdss3.org/collaboration/boiler-plate.php |
| 2MASS | 37 | in the Note in full |

A version with the DESI and short Legacy Surveys texts in full and the body cut
(`~/star/paper/rnaas/debbie_rnaas_fullack`) is 1,526 words, 26 over, and still leaves out DES
and SDSS-III. Not for submission.
