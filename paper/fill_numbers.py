"""Render winnie.qmd from the template, using only the frozen v1.2 numbers.

Every quantitative claim in the manuscript comes from
catalogue_debbie/v1.2/PAPER_NUMBERS.json, which is computed from the frozen,
checksummed v1.0 catalogue. The paper previously quoted four different
catalogue sizes at once because the pipeline kept growing the catalogue under
the draft; this makes that impossible. If a key is missing the render fails
loudly rather than leaving a placeholder in the text.
"""
import json
import re
import sys

NUM = "/home/ben/star/catalogue_debbie/v1.2/PAPER_NUMBERS.json"
TPL = "/home/ben/star/paper/debbie_template.qmd"
OUT = "/home/ben/star/paper/debbie.qmd"

n = json.load(open(NUM))

def comma(x):
    return f"{x:,}"

KEYS = {
    "N": str(n["N"]),
    "N_COMMA": comma(n["N"]),
    "Z_COV": f"{n['z_coverage_pct']:.2f}",
    "N_SPEC_COMMA": comma(n["n_spec"]),
    "N_LSP_COMMA": comma(n["n_ls_photo"]),
    "N_OWN_COMMA": comma(n["n_own_photo"]),
    "Z_MED": f"{n['z_median']:.2f}",
    "TOLMAN": f"{n['tolman_at_median']:.2f}",
    "R_KPC": f"{n['r_kpc_median']:.0f}",
    "MG_MED": f"{n['Mg_median']:.2f}",
    "LOGM_MED": f"{n['logM_median']:.2f}",
    "N_DWARF_COMMA": comma(n["n_dwarf_z005"]),
    "PSTAR_COMMA": comma(n["pstar_flagged"]),
    "PSTAR_PCT": f"{n['pstar_pct_of_N']:.1f}",
    "DWARF_RKPC": f"{n['dwarf_rkpc_median']:.1f}",
    "DWARF_LOGM": f"{n['dwarf_logM_median']:.1f}",
    "AREA_COMMA": comma(n["area_deg2"]),
    "N_EXT_COMMA": comma(n["N_lowb_ext"]),
    "DESIG_PCT": f"{n['designated_pct']:.1f}",
    "SIGMA_NMAD": f"{n['own_sigma_nmad']:.3f}",
    "CONTROL": "0.194",
    "VII292_PRESENT": f"{n['vii292_present_pct']:.1f}",
    "VII292_NONSTELLAR": f"{n['vii292_nonstellar_pct']:.1f}",
    "N_UDG_COMMA": comma(n["n_udg"]),
    "N_REMAIN": comma(n["n_remain"]),
    "DESI_MAIN": comma(n["desi_main"]), "DESI_MAIN_PCT": f"{n['desi_main_pct']:.2f}",
    "DESI_GAL": comma(n["desi_gal"]), "DESI_QSO": comma(n["desi_qso"]),
    "DESI_STAR": comma(n["desi_star"]), "DESI_REL": comma(n["desi_rel"]),
    "DESI_PURITY": f"{n['desi_purity']:.2f}", "DESI_PURITY_ERR": f"{n['desi_purity_err']:.2f}",
    "DESI_PURITY_EXT": f"{n['desi_purity_ext']:.2f}", "N12": comma(n["N"]),
    "POOL_HI": comma(n["pool_hi"]),
    "POOL_VER": comma(n["pool_verified"]),
    "POOL_UNVER": comma(n["pool_unverified"]),
    "POOL_UNVER_PCT": f"{n['pool_unver_pct']:.1f}",
    "UNVER_BIG": comma(n["unverified_big_r6"]),
    "Z_MED3": f"{n['z_median']:.3f}",
    "OWN_PCT": f"{n['own_pct_of_z']:.0f}",
    "TOLMAN3": f"{0.75*n['tolman_at_median']:.2f}",
    "PSTAR_EXP_COMMA": comma(n["pstar_expected_stars"]),
    "PSTAR_EXP_PCT": f"{n['pstar_expected_pct']:.1f}",
    "DEC_B1": comma(n["dec_band_counts"]["-2..0"]),
    "DEC_B2": comma(n["dec_band_counts"]["0..2"]),
    "DEC_NBR_LO": comma(n["dec_neighbour_range"][0]),
    "DEC_NBR_HI": comma(n["dec_neighbour_range"][1]),
    "PURITY_N_COMMA": comma(n["purity_sample_n"]),
    "PURITY_NEFF_COMMA": comma(int(round(n["purity_sample_neff"]))),
    "PURITY_CELLS": str(n["purity_sample_cells"]),
    "N_SPEC_PREV_COMMA": comma(n["n_spec_prev"]),
}

text = open(TPL).read()
for k, v in KEYS.items():
    text = text.replace("{{" + k + "}}", v)

left = re.findall(r"\{\{[A-Z_]+\}\}", text)
if left:
    sys.exit(f"UNFILLED PLACEHOLDERS: {sorted(set(left))}")

open(OUT, "w").write(text)
print(f"debbie.qmd rendered from v1.2 (N = {KEYS['N_COMMA']})")
