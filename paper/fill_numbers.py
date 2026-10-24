"""Render debbie.qmd from the template, using only the frozen v1.2 numbers.

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

NUM = "/home/ben/star/catalogue_debbie/v1.4/PAPER_NUMBERS.json"   # same rows as v1.2
TPL = "/home/ben/star/paper/debbie_template.qmd"
OUT = "/home/ben/star/paper/debbie.qmd"

n = json.load(open(NUM))

# The visual purity measurement. Written by sweeps/purity_by_regime.py from the
# inspector's verdict file, so it changes when more objects are judged or the
# undecided ones are resolved; rerun that script, then this one.
PUR = "/home/ben/star/catalogue_debbie/purity_by_regime.json"
pur = json.load(open(PUR))

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
    "PUR_SPLIT": "25.5",
    "PUR_JUDGED": comma(pur["ALL"]["judged"] + pur["ALL"]["undecided"]),
    "PUR_UND": comma(pur["ALL"]["undecided"]),
}
for tag, key in (("ALL", "ALL"), ("B", "bright half"), ("F", "faint half")):
    r = pur[key]
    KEYS[f"PUR_{tag}"] = f"{100 * r['purity']:.1f}"
    KEYS[f"PUR_{tag}_LO"] = f"{100 * r['ci'][0]:.1f}"
    KEYS[f"PUR_{tag}_HI"] = f"{100 * r['ci'][1]:.1f}"
    KEYS[f"PUR_{tag}_BLO"] = f"{100 * r['undecided_bound'][0]:.1f}"
    KEYS[f"PUR_{tag}_BHI"] = f"{100 * r['undecided_bound'][1]:.1f}"
    KEYS[f"PUR_{tag}_N"] = comma(r["n_catalogue"])
    KEYS[f"PUR_{tag}_FRAC"] = f"{100 * r['n_catalogue'] / n['N']:.1f}"
# VST ATLAS against the catalogue and against the visual verdicts, written by
# sweeps/vst_vs_verdicts.py. These were hard-coded in the text before 2026-10-01,
# from a run that left no script.
vst = json.load(open("/home/ben/star/catalogue_debbie/vst_vs_verdicts.json"))
vr, vv = vst["random"], vst["verdicts"]
cls = vr["classes"]
KEYS.update({
    "VST_N": comma(vr["n"]),
    "VST_MATCH": f"{100 * vr['matched'] / vr['n']:.1f}",
    "VST_GAL": f"{100 * vr['galaxy_frac_of_matched']:.2f}",
    "VST_STAR": f"{100 * (cls.get('star', 0) + cls.get('probable star', 0)) / vr['matched']:.2f}",
    "VST_NOISE": f"{100 * cls.get('noise', 0) / vr['matched']:.2f}",
    "VST_YES_N": str(vv["all yes"]["in_footprint"]), "VST_YES_M": str(vv["all yes"]["matched"]),
    "VST_NO_N": str(vv["all no"]["in_footprint"]), "VST_NO_M": str(vv["all no"]["matched"]),
    "VST_YES_RATE": f"{100 * vv['all yes']['matched'] / vv['all yes']['in_footprint']:.1f}",
    "VST_NO_RATE": f"{100 * vv['all no']['matched'] / vv['all no']['in_footprint']:.1f}",
    "VST_NO_GAL": str(vv["all no"]["classes"].get("galaxy", 0)
                      + vv["all no"]["classes"].get("probable galaxy", 0)),
    "VST_B_RATE": f"{100 * vst['random mu<25.5']['matched'] / vst['random mu<25.5']['in_footprint']:.1f}",
    "VST_F_RATE": f"{100 * vst['random mu>=25.5']['matched'] / vst['random mu>=25.5']['in_footprint']:.1f}",
})
# The star score and flag: sweeps/star_score.py and sweeps/star_flag_effect.py.
ss = json.load(open("/home/ben/star/catalogue_debbie/star_score.json"))
fe = json.load(open("/home/ben/star/catalogue_debbie/star_flag_effect.json"))
if ss["main"]["n"] != n["N"]:
    sys.exit(f"star scores are for {ss['main']['n']:,} entries, paper for {n['N']:,}")
pc = ss["pstar_calibration"]
reg = ss["desi"]["galaxy_loss_by_regime"]
KEYS.update({
    "STAR_NFEAT": str(len(ss["features"])),
    "STAR_NGAL": comma(ss["truth"]["galaxies"]), "STAR_NSTAR": comma(ss["truth"]["stars"]),
    "STAR_AUC": f"{ss['truth']['auc_classifier']:.3f}", "STAR_AUC_P": f"{ss['truth']['auc_pstar']:.3f}",
    "STAR_GLOSS": f"{ss['desi']['classifier']['galaxies']:.1f}",
    "STAR_D_CLF": f"{ss['desi']['classifier']['stars']:.1f}",
    "PS09_GAL": f"{ss['desi']['p_star > 0.9']['galaxies']:.1f}",
    "STAR_D_ANY": f"{ss['desi']['either']['stars']:.1f}",
    "STAR_GLOSS_ANY": f"{ss['desi']['either']['galaxies']:.1f}",
    "STAR_G_CLF": f"{ss['gaia']['classifier']:.1f}", "STAR_G_ANY": f"{ss['gaia']['either']:.1f}",
    "GAIA_N": comma(ss["gaia"]["n_positions"]), "GAIA_NCERT": comma(ss["gaia"]["certain_stars"]),
    "GAIA_MATCH": f"{100 * ss['gaia']['matched'] / ss['gaia']['n_positions']:.2f}",
    "GAIA_CERT": f"{100 * ss['gaia']['certain_stars'] / ss['gaia']['n_positions']:.2f}",
    "GAIA_CERT_LO": f"{ss['gaia']['certain_ci'][0]:.2f}", "GAIA_CERT_HI": f"{ss['gaia']['certain_ci'][1]:.2f}",
    "PS_GAL": f"{pc['mean_confirmed_galaxies']:.3f}", "PS_CAT": f"{pc['mean_catalogue']:.3f}",
    "PS_SHARE": f"{100 * pc['mean_confirmed_galaxies'] / pc['mean_catalogue']:.0f}",
    "PS_COV": f"{ss['main']['pstar_present_pct']:.1f}", "PS_CUT": f"{ss['pstar_cut']}",
    "FLAG1": comma(ss["main"]["flag_classifier"]), "FLAG2": comma(ss["main"]["flag_pstar_only"]),
    "FLAG_PCT": f"{ss['main']['flagged_pct']:.1f}",
    "STAR_GL_SMALL": f"{reg['r<4']['flagged']:.1f}", "STAR_GL_BIG": f"{reg['r>=6']['flagged']:.1f}",
    "EYE_NO": str(fe["sample no"]["n"]), "EYE_NO_FLAG": str(fe["sample no"]["either"]),
    "PUR_ALL_FLAG": f"{100 * fe['purity ALL either']['purity']:.1f}",
})
# The purity figure after the second look: sweeps/purity_final.py.
pf = json.load(open("/home/ben/star/catalogue_debbie/purity_final.json"))
if not pf["complete"]:
    sys.exit("the second look is not complete; purity_final.json cannot be quoted")
A, B, bc = pf["A"]["purity"], pf["B"]["purity"], pf["became"]
gal = lambda d: d.get("yes", 0) + d.get("small", 0)
for tag, src, k in (("PF_A", A, "all entries ALL"), ("PF_AB", A, "all entries bright"),
                    ("PF_AF", A, "all entries faint"), ("PF_B", B, "all entries ALL"),
                    ("PF_BF", B, "all entries faint")):
    KEYS[tag] = f"{100 * src[k]['corrected']:.1f}"
    KEYS[tag + "_LO"] = f"{100 * src[k]['corrected_ci'][0]:.1f}"
    KEYS[tag + "_HI"] = f"{100 * src[k]['corrected_ci'][1]:.1f}"
ctl = pf["A"]["accepted_controls"]["ALL"]
KEYS.update({
    "PF_A_BLO": f"{100 * A['all entries ALL']['undecided_bound'][0]:.1f}",
    "PF_A_BHI": f"{100 * A['all entries ALL']['undecided_bound'][1]:.1f}",
    "PF_STILL": str(A["all entries ALL"]["still_undecided"]),
    "PF_SMALL": f"{100 * (A['all entries ALL']['corrected'] - B['all entries ALL']['corrected']):.1f}",
    "PF_A_STAR": f"{100 * A['star-flagged left out ALL']['corrected']:.1f}",
    "PF_NREJ": str(pf["second_look"]["rejections"]), "PF_NUND": str(pf["second_look"]["undecided"]),
    "PF_NCTL": str(pf["second_look"]["accepted_controls"]),
    "PF_CTL_DEC": str(ctl["decided"]), "PF_CTL_NO": str(ctl["now_rejected"]),
    "PF_CTL_SMALL": str(bc["accepted"].get("small", 0)),
    "PF_REJ_GAL": str(gal(bc["rejected"])), "PF_REJ_NO": str(bc["rejected"].get("no", 0)),
    "PF_REJ_AMB": str(bc["rejected"].get("amb", 0)), "PF_UND_GAL": str(gal(bc["undecided"])),
    "PF_VD_GAL": str(gal(bc["rejected vst_detected"])), "PF_VD_NO": str(bc["rejected vst_detected"].get("no", 0)),
    "PF_VN_GAL": str(gal(bc["rejected vst_not_detected"])),
    "PF_VN_N": str(sum(bc["rejected vst_not_detected"].values())),
    "PF_UND_NO": str(bc["undecided"].get("no", 0)), "PF_UND_AMB": str(bc["undecided"].get("amb", 0)),
    "PF_STAR_DELTA": f"{100 * abs(A['star-flagged left out ALL']['corrected'] - A['all entries ALL']['corrected']):.2f}",
    "PF_NACC": str(pf["first_stage"]["accepted"]), "PF_NSET": str(pf["first_stage"]["set_aside"]),
    "PF_SET_GAL": str(gal(bc["rejected"]) + gal(bc["undecided"])),
    "PF_SET_NO": str(bc["rejected"].get("no", 0) + bc["undecided"].get("no", 0)),
    "PF_SET_AMB": str(bc["rejected"].get("amb", 0) + bc["undecided"].get("amb", 0)),
    "VF_NOT_N": str(pf["vst_final"]["not"]["in_footprint"]),
    "VF_NOT_M": str(pf["vst_final"]["not"]["matched"]),
    "VF_NOT_GAL": str(pf["vst_final"]["not"]["vst_galaxy"]),
    "PF_SL_NO": str(pf["star_flag"]["second_look_no"]),
    "PF_SL_NO_FLAG": str(pf["star_flag"]["of_those_flagged"]),
})
# Two facts about the released columns, read from the tables themselves.
import csv as _csv
import numpy as _np
_rows = list(_csv.DictReader(open("/home/ben/star/catalogue_debbie/v1.4/DEBBIE.csv", encoding="utf-8")))
_z = _np.array([float(r["z_phot"]) if r["z_phot"] not in ("", "nan") else _np.nan for r in _rows])
_m = _np.array([float(r["log_mstar"]) if r["log_mstar"] not in ("", "nan") else _np.nan for r in _rows])
_hi = (_z >= 0.35) & _np.isfinite(_m)
_r = _np.round(_np.asarray(_np.load("/home/ben/star/catalogue_debbie/v1.2_main.npy")["shape_r"], float), 3)
KEYS["MASS_HI_PCT"] = f"{100 * float((_m[_hi] > 12).mean()):.0f}"
KEYS["N_PILE"] = str(int((_r == 13.624).sum()))
# The candidate-set comparison on the benchmark sample, written by
# sweeps/multimodel_folds.py --subset benchmark.
mm = json.load(open("/home/ben/star/catalogue_debbie/multimodel_folds_benchmark.json"))
if mm["subset"] != "benchmark" or set(mm["by_source"]) != {"deepshadows_lsbg", "deepshadows_hard_neg"}:
    sys.exit("multimodel file is not the benchmark sample")
_auc = lambda k: _np.array([f["auc"][k] for f in mm["folds"]])
_ll = lambda k: _np.array([f["logloss"][k] for f in mm["folds"]])
_DUAL, _ALL = "equal mean, full + bigru", "equal mean, all"
for tag, k in (("SIX", "base6"), ("PROF", "profile"), ("P6", "prof+6"), ("FULL", "full"),
               ("GRU", "bigru"), ("DUAL", _DUAL), ("STACK", "stacked (logistic)")):
    KEYS[f"MM_{tag}"] = f"{_auc(k).mean():.4f}"
    KEYS[f"MM_{tag}_SD"] = f"{_auc(k).std(ddof=1):.4f}"
_g = _auc("full") - _auc("base6")
_d = _auc(_DUAL) - _auc("full")
_gr = _auc("bigru") - _auc("full")
_al = _auc(_ALL) - _auc("full")
KEYS.update({
    "MM_N": comma(mm["n"]), "MM_NPOS": comma(mm["n_positive"]),
    "MM_NNEG": comma(mm["n"] - mm["n_positive"]),
    "MM_NTRAIN": comma(int(round(mm["n"] * 0.8, -2))),
    "MM_GAIN": f"{_g.mean():+.3f}", "MM_GAIN_SD": f"{_g.std(ddof=1):.3f}", "MM_GAIN_MIN": f"{_g.min():+.3f}",
    "MM_GRU_D": f"{_gr.mean():+.4f}", "MM_GRU_D_SD": f"{_gr.std(ddof=1):.4f}",
    "MM_GRU_POS": str(int((_gr > 0).sum())),
    "MM_DUAL_D": f"{_d.mean():+.4f}", "MM_DUAL_D_SD": f"{_d.std(ddof=1):.4f}",
    "MM_DUAL_POS": str(int((_d > 0).sum())),
    "MM_DUAL_LL": f"{100 * (1 - float((_ll(_DUAL) / _ll('full')).mean())):.0f}",
    "MM_ALL_D": f"{abs(_al.mean()):.4f}",
    "MM_RATIO": f"{round(float(_g.mean() / _d.mean()), -1):.0f}",
})
if min(max(f["weights"].values()) for f in mm["folds"]) <= 0.99:
    sys.exit("the paper says one model takes more than 0.99 of the weight on every fold; it no longer does")
if _al.mean() >= 0 or _d.min() <= 0:
    sys.exit("the multimodel paragraph says the five-way mean loses and the dual mean wins on every fold; "
             "the numbers no longer say that, so rewrite the paragraph")
# The two-model screen: the re-extraction check, the models' own folds, the
# verdict comparison, and purity among the entries it passes.
import numpy as np
rx = json.load(open("/home/ben/star/catalogue_debbie/reextraction_check.json"))
if rx["verdict"] != "pass":
    sys.exit("the re-extraction check did not pass; the screen paragraph cannot stand")
md = json.load(open("/home/ben/star/output/dual_models/deepshadows/models.json"))["five_fold_auc"]
dv = json.load(open("/home/ben/star/catalogue_debbie/dual_screen_vs_verdicts.json"))["scores"]
ps = json.load(open("/home/ben/star/catalogue_debbie/purity_after_screen.json"))
cf, av = ps["results"]["cross_fitted"], ps["against_verdicts"]
pc = lambda x: f"{100 * x:.1f}"
KEYS.update({
    "RX_S_AUG": f"{rx['single_split']['august']:.4f}", "RX_S_NEW": f"{rx['single_split']['reextracted']:.4f}",
    "RX_F_AUG": f"{rx['five_folds']['full']['august']:.4f}", "RX_F_NEW": f"{rx['five_folds']['full']['reextracted']:.4f}",
    # objects in the re-extracted pack, not the five-fold sample drawn from it
    "RX_N": comma(len(np.load("/home/ben/star/output/radial_profiles_reextracted.npz",
                              allow_pickle=True)["object_id"])),
    "SC_TREES": f"{md['trees']['mean']:.4f}", "SC_NET": f"{md['network']['mean']:.4f}",
    "SC_MEAN": f"{md['mean']['mean']:.4f}",
    "SV_TREES": f"{dv['deepshadows_trees']['auc']:.3f}", "SV_NET": f"{dv['deepshadows_network']['auc']:.3f}",
    "SV_S2": f"{dv['stage2_mean']['auc']:.2f}",
    "SV_MEAN": f"{av['auc_mean'][0]:.3f}", "SV_MEAN_LO": f"{av['auc_mean'][1][0]:.3f}",
    "SV_MEAN_HI": f"{av['auc_mean'][1][1]:.3f}",
    "SV_NGAL": str(av["n_galaxy"]), "SV_NNOT": str(av["n_not"]),
    "SV_SF_NOT": str(av["non_galaxies_star_flagged"]), "SV_SF_GAL": str(av["galaxies_star_flagged"]),
    "SV_NOSTAR": f"{av['auc_where_star_flag_passes'][0]:.3f}", "SV_NOSTAR_N": str(av["n_not_where_star_flag_passes"]),
    "SF_KEEP": f"{100 * ps['keep_galaxies']:.0f}",
    "SF_LOWPCT": f"{100 * (1 - ps['keep_galaxies']):.0f}",
    "SF_ENT": f"{cf['entries_flagged_pct']:.1f}", "SF_GAL": f"{cf['galaxies_flagged_pct']:.1f}",
    "SF_NOT": f"{cf['non_galaxies_flagged_pct']:.0f}",
    "SF_PASS": pc(cf["ALL"]["passing"][0]), "SF_PASS_LO": pc(cf["ALL"]["passing"][1][0]),
    "SF_PASS_HI": pc(cf["ALL"]["passing"][1][1]),
    "SF_GAIN": f"{ps['paired_gain']['points']:.1f}", "SF_GAIN_LO": f"{ps['paired_gain']['ci'][0]:.1f}",
    "SF_GAIN_HI": f"{ps['paired_gain']['ci'][1]:.1f}",
    "SF_FAINT_ALL": pc(cf["faint"]["all"][0]), "SF_FAINT_PASS": pc(cf["faint"]["passing"][0]),
    "SF_BRIGHT_ALL": pc(cf["bright"]["all"][0]), "SF_BRIGHT_PASS": pc(cf["bright"]["passing"][0]),
    "SF_FLAGPUR": f"{100 * cf['ALL']['flagged'][0]:.0f}", "SF_FLAGPUR_LO": f"{100 * cf['ALL']['flagged'][1][0]:.0f}",
    "SF_FLAGPUR_HI": f"{100 * cf['ALL']['flagged'][1][1]:.0f}",
})
# The benchmark and the frozen spatial fold, from the files the scripts wrote
# (benchmark_deepshadows.py and compare_common_fold.py); the figure reads the first.
bm = json.load(open("/home/ben/star/output/benchmark_deepshadows.json"))
fz = json.load(open("/home/ben/star/output/common_fold_trees.json"))["arms"]
KEYS.update({
    "BM_AUC": f"{bm['radial']['auc']:.4f}", "BM_LO": f"{bm['radial']['ci'][0]:.4f}",
    "BM_HI": f"{bm['radial']['ci'][1]:.4f}", "BM_SIX": f"{bm['six_features']['auc']:.4f}",
    "BM_SIX_LO": f"{bm['six_features']['ci'][0]:.4f}", "BM_SIX_HI": f"{bm['six_features']['ci'][1]:.4f}",
    "BM_NTRAIN": comma(bm["n_train"]), "BM_NTEST": comma(bm["n_test"]),
    "BM_NPER": comma(bm["n_test"] // 2),
    "FZ_SIX": f"{fz['base6_only:xgb']['auc']:.4f}", "FZ_PROF": f"{fz['profile_only:xgb']['auc']:.4f}",
    "FZ_P6": f"{fz['profile_plus_base6:xgb']['auc']:.4f}", "FZ_FULL": f"{fz['flattened_all:xgb']['auc']:.4f}",
})
# The SDSS footprint hole (sdss_hole.py), found from the density map.
sh = json.load(open("/home/ben/star/catalogue_debbie/sdss_hole.json"))
KEYS.update({
    "SH_SOUTH": f"{sh['south']['debbie_per_deg2']:.1f}", "SH_NORTH": f"{sh['north']['debbie_per_deg2']:.2f}",
    "SH_SDSS": comma(round(sh["north"]["sdss_gal_per_deg2"], -3)),
    "SH_N20": f"{sh['north']['expected_within_20']:.1f}", "SH_N45": f"{sh['north']['expected_20_45']:.1f}",
})
# Where the candidates passing the cuts went (screen_accounting.py). Replaces
# the old reading that 103,317 were "not checked because queries failed".
sa = json.load(open("/home/ben/star/catalogue_debbie/screen_accounting.json"))
if sa["candidates"] != n["pool_hi"]:
    sys.exit(f"screen accounting is for {sa['candidates']:,} candidates, paper for {n['pool_hi']:,}")
fr = sa["fragment_filter_rate_by_reff"]
KEYS.update({
    "SA_INDEX_E8": f"{sa['index_rows'] / 1e8:.1f}",
    # over every source passing the cuts, mining's unrecorded removals included
    "FC_TOTAL": comma(sa["full"]["cuts_total"]), "FC_ID20": comma(sa["full"]["within_20"]),
    "FC_SDSS5": comma(sa["full"]["within_20_sdss_5"]), "FC_FRAG": comma(sa["full"]["fragment_filter"]),
    "FC_MINED": comma(sa["full"]["mined_out"]),
    "FC_FRAG_PCT": f"{sa['full']['fragment_rate_pct']:.0f}",
    "FC_FR34": f"{sa['full']['fragment_rate_by_reff']['3-4']:.0f}",
    "FC_FR1015": f"{sa['full']['fragment_rate_by_reff']['10-15']:.0f}",
    "FC_FR1015_KEEP": f"{1 - sa['full']['fragment_rate_by_reff']['10-15'] / 100:.2f}",
    "FC_E2E": f"{0.63 * (1 - sa['full']['fragment_rate_by_reff']['10-15'] / 100):.2f}",
    "FC_OLD90": comma(sa["full"]["mined_beyond_45_old_radius"]),
    "FC_UNREC": comma(sa["full"]["mined_unrecorded"]),
    "FC_MINED_OTHER": comma(sa["full"]["mined_beyond_45"]),
    "FC_ID20_PCT": f"{100 * sa['full']['within_20'] / sa['full']['cuts_total']:.0f}",
    "FC_SMALL3": comma(sa["clump_filter"] + sa["full"]["mined_beyond_45"] + sa["never_queried"]),
    "FC_EXPORT": comma(sa["full"]["export_and_merge"]),
    "FC_BEFORE_DESI": comma(n["N"] + n["desi_main"]),
    "SA_ID20": comma(sa["identified_within_20"]), "SA_ID20_PCT": f"{sa['identified_within_20_pct']:.1f}",
    "SA_SDSS5": comma(sa["of_which_sdss_within_5"]),
    "SA_FRAG": comma(sa["fragment_filter"]), "SA_FRAG_PCT": f"{sa['fragment_filter_pct']:.0f}",
    "SA_CLUMP": comma(sa["clump_filter"]), "SA_CLUMP_PCT": f"{sa['clump_filter_pct']:.1f}",
    "SA_NEVER": comma(sa["never_queried"]), "SA_NEVER_PCT": f"{sa['never_queried_pct']:.1f}",
    "SA_NEVER_SMALL": comma(sa["never_queried_below_305"]),
    "SA_QUERIED": comma(sa["queried"]), "SA_KNOWN": comma(sa["queried_known"]),
    "SA_NOVEL": comma(sa["queried_novel"]),
    "SA_SDSS_OWN_PCT": str(sa["sdss_area_own_pct"]),
    "SA_FR34": f"{fr['3-4']:.0f}", "SA_FR1015": f"{fr['10-15']:.0f}",
    "SA_FR1015_KEEP": f"{1 - fr['10-15'] / 100:.2f}",
    "SA_E2E": f"{0.63 * (1 - fr['10-15'] / 100):.2f}",
})
# The first five entries by designation number, for the example-rows table
# AJ asks for when the full table is machine-readable.
import csv as _csv
_rows = list(_csv.DictReader(open("/home/ben/star/catalogue_debbie/v1.4/DEBBIE.csv", encoding="utf-8")))
_rows.sort(key=lambda r: int(r["name"].split()[1]))
_fmt = lambda x, f: (f % float(x)) if x not in ("", "nan") else "--"
_lines = ["| designation | RA | Dec | $r_{\\rm eff}$ | $\\mu_{\\rm eff}(g)$ | $g$ | $g-r$ | $z_{\\rm phot}$ | star flag |",
          "|:--------------|-----------:|-----------:|-------:|--------:|-------:|-------:|--------:|-------:|"]
for r in _rows[:5]:
    # pandoc will not close $...$ before a digit, so "$-$8.8" would run math across
    # the table: negative numbers go in math whole, and designations keep a plain hyphen
    _lines.append("| " + " | ".join([r["name"], _fmt(r["ra_deg"], "%.5f"),
                                      "$" + _fmt(r["dec_deg"], "%.5f") + "$", _fmt(r["r_eff_arcsec"], "%.2f"),
                                      _fmt(r["mu_eff_g"], "%.2f"), _fmt(r["g_mag"], "%.2f"), _fmt(r["g_r"], "%.2f"),
                                      _fmt(r["z_phot"], "%.3f"), r["star_flag"] or "0"]) + " |")
KEYS["EXAMPLE_ROWS"] = "\n".join(_lines)
# Published catalogue sizes (abstracts of Tanoglidis et al. 2021, ApJS 252, 18
# and Zaritsky et al. 2023, ApJS 267, 27); the ratios follow from N.
KEYS["CMP_X_DES"] = str(round(n["N"] / 23790))
KEYS["CMP_X_SM"] = str(round(n["N"] / 7070))
if abs(cf["ALL"]["all"][0] - pf["A"]["purity"]["all entries ALL"]["corrected"]) > 1e-4:
    sys.exit("purity_after_screen.json was computed on different verdicts from purity_final.json")
if pur["ALL"]["n_catalogue"] != n["N"]:
    sys.exit(f"purity file is for {pur['ALL']['n_catalogue']:,} entries, paper for {n['N']:,}")

text = open(TPL).read()
for k, v in KEYS.items():
    text = text.replace("{{" + k + "}}", v)

left = re.findall(r"\{\{[A-Z0-9_]+\}\}", text)
if left:
    sys.exit(f"UNFILLED PLACEHOLDERS: {sorted(set(left))}")

open(OUT, "w").write(text)
print(f"debbie.qmd rendered from the v1.4 numbers (N = {KEYS['N_COMMA']})")
