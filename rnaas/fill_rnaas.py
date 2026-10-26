"""Render the DEBBIE Research Note from its template, every number from a file.

Ben, 2026-10-20: "let's do RNAAS first", then the full paper in AJ. C. Lintott
(RNAAS, 2026-10-20) allows a Note before the journal paper, and asks that the
Note give enough of the selection for the catalogue to be usable before the
paper exists.

Sources:
  counts, numbers, parts     catalogue_debbie/v1.5/numbers.json, numbers.npz
  index size                 screen_accounting.json, wide_screen_count.json
  area, median photo-z       v1.4/PAPER_NUMBERS.json (the original part; n_spec = 0,
                             so every redshift is Duncan 2022 or own_photoz.py)
  selection function         release/catalogue/selection_function.json
  purity                     purity_final.json (original), purity_v3.json
                             "ruled + rule" (recovered = add, extension = ext)
  scores                     output/benchmark_deepshadows.json,
                             dual_screen_vs_verdicts.json (first sample),
                             v3_model_check.json (second sample, v1.5 members)
  purity by neighbour, score purity_by_property.json
  astrometry, scored         v1.5/cds_entry.json (sweeps/cds_entry_v15.py)

A key missing from KEYS, or one left in the output, stops the render. The
word count (texcount, references included) and the abstract's count are
printed against the AAS limits of 1,500 and 150.

    .venv/bin/python paper/rnaas/fill_rnaas.py [--fullack]
"""
import json
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# --fullack: the variant with the DESI DR1 and Legacy Surveys (short)
# acknowledgements in full and the body cut to fit, in case the AAS counts them
FULLACK = "--fullack" in sys.argv
TPL = f"{HERE}/debbie_rnaas_template{'_fullack' if FULLACK else ''}.tex"
OUT = f"{HERE}/debbie_rnaas{'_fullack' if FULLACK else ''}.tex"
C = "/home/ben/star/catalogue_debbie"
TEXCOUNT = "/home/ben/.TinyTeX/bin/x86_64-linux/texcount"
ZENODO = r"\textbf{DOI to come}"          # after CDS enters DEBBIE in the Dictionary


def load(p):
    return json.load(open(p if p.startswith("/") else f"{C}/{p}"))


def comma(x):
    return f"{x:,}"


num = load("v1.5/numbers.json")
groups = {g["group"]: g for g in num["groups"]}
pn = load("v1.4/PAPER_NUMBERS.json")
if pn["n_spec"] != 0:                    # the abstract says every redshift is photometric
    raise SystemExit(f"{pn['n_spec']} spectroscopic redshifts; the abstract says none")
sf = load("/home/ben/star/release/catalogue/selection_function.json")
pf = load("purity_final.json")
p3 = load("purity_v3.json")["variants"]["ruled + rule"]
bm = load("/home/ben/star/output/benchmark_deepshadows.json")
dv = load("dual_screen_vs_verdicts.json")["scores"]
v3 = load("v3_model_check.json")["results"]
ce = load("v1.5/cds_entry.json")
pp = load("purity_by_property.json")             # sweeps/purity_by_property.py
sample3 = load("/home/ben/legacy/inspect_sample_v3.json")

# Selection function: cells inside the window, 3-15", with at least 30 known
# galaxies (as the full paper quotes it); and every window cell above 15".
WIN = ("24-25", "25-26", "26-26.5")
SMALL = ("3-4", "4-5", "5-6", "6-8", "8-10", "10-15")
cells = [(k.split("_")[0], k.split("_")[1], v) for k, v in sf["grid"].items()]
small = [v["frac"] for r, m, v in cells if r in SMALL and m in WIN and v["n"] >= 30]
big = [v for r, m, v in cells if r == "15-100" and m in WIN]
big_n = sum(v["n"] for v in big)
big_k = sum(round(v["frac"] * v["n"]) for v in big)

# Purity by part: (left out, its interval, undecided counted as no)
o = pf["A"]["purity"]["all entries ALL"]
pur = {"original": (o["corrected"], o["corrected_ci"], o["undecided_bound"][0]),
       "recovered": (p3["add"]["left_out"], p3["add"]["left_out_ci"], p3["add"]["as_no"]),
       "extension": (p3["ext"]["left_out"], p3["ext"]["left_out_ci"], p3["ext"]["as_no"])}
auc = {"original": dv["deepshadows_mean"]["auc"],
       "recovered": v3["add screened"]["deepshadows_mean"][0],
       "extension": v3["ext screened"]["deepshadows_mean"][0]}

# Density inside and just south of the SDSS footprint, 150 < RA < 200.
z = np.load(f"{C}/v1.5/numbers.npz", allow_pickle=True)


def density(d0, d1):
    m = (z["ra"] > 150) & (z["ra"] < 200) & (z["dec"] > d0) & (z["dec"] < d1)
    area = np.degrees(np.radians(50)) * np.degrees(np.sin(np.radians(d1)) - np.sin(np.radians(d0)))
    return m.sum() / area


pct = lambda x: f"{100 * x:.1f}"                                         # noqa: E731
BAND = {"original": r"$\ge 30^\circ$", "recovered": r"$\ge 30^\circ$", "extension": r"12--30$^\circ$"}
rows = []
for part in ("original", "recovered", "extension"):
    p, ci, as_no = pur[part]
    rows.append(f"{part.capitalize()} & {BAND[part]} & {comma(num['by_sub'][part])} & "
                f"{pct(p)} ({pct(ci[0])}--{pct(ci[1])}) & {pct(as_no)} & {auc[part]:.2f} & "
                f"{pct(ce['by_sub'][part]['within_1'])} \\\\")
rows.append(f"All & & {comma(num['entries'])} & & & & {pct(ce['astrometry']['within_1'])} \\\\")

# Where purity is high and low (second sample by distance to a catalogued
# galaxy; every part above a mean score of 0.2). The Note says the fixed radii
# kept catalogued galaxies out of the original part: check it on its sample.
near = pp["second"]["nearest_galaxy_split"]
if pp["first_min_galaxy_sep"] < 45:
    raise SystemExit(f"an original sample object lies {pp['first_min_galaxy_sep']}\" from a catalogued galaxy")
above = {part: pp["score"][f"{part} >= 0.2"] for part in ("original", "recovered", "extension")}
others = (above["recovered"]["share_of_part_above"] * num["by_sub"]["recovered"]
          + above["extension"]["share_of_part_above"] * num["by_sub"]["extension"]) \
    / (num["by_sub"]["recovered"] + num["by_sub"]["extension"])
lo02 = min(a["above"]["purity"] for a in above.values())
hi02 = max(a["above"]["purity"] for a in above.values())

n_gal = load("wide_screen_count.json")["index_objects_used"]
KEYS = {
    "N": comma(num["entries"]),
    "INDEX_GAL_E8": f"{n_gal / 1e8:.1f}",
    "INDEX_GAL": comma(n_gal),
    "THR": f"{dv['deepshadows_mean']['threshold_keeping_95pct_of_galaxies']:.2f}",
    "THR_NOT": f"{dv['deepshadows_mean']['non_galaxies_below_it_pct']:.0f}",
    "N_DWARF": comma(pn["n_dwarf_z005"]),
    "DWARF_PCT": f"{100 * pn['n_dwarf_z005'] / pn['N']:.1f}",
    "AREA": comma(pn["area_deg2"]),
    "N_EXT1": comma(groups["ext_v1.2"]["kept"]),
    "PUR_LO": f"{100 * min(v[0] for v in pur.values()):.0f}",
    "PUR_HI": f"{100 * max(v[0] for v in pur.values()):.0f}",
    "PUR_O": pct(pur["original"][0]), "PUR_R": pct(pur["recovered"][0]), "PUR_E": pct(pur["extension"][0]),
    "PUR1_N": comma(pf["n_first"]),
    "PUR3_N": comma(load("purity_v3.json")["n_judged"]),
    "PUR3_DRAWN": comma(len(sample3["objects"])),
    "SF_LO": f"{100 * min(small):.0f}", "SF_HI": f"{100 * max(small):.0f}",
    "SF_BIG_K": str(big_k), "SF_BIG_N": str(big_n),
    "BM_AUC": f"{bm['radial']['auc']:.3f}",
    "AUC_O": f"{auc['original']:.2f}", "AUC_R": f"{auc['recovered']:.2f}", "AUC_E": f"{auc['extension']:.2f}",
    "AUC_S2_O": f"{dv['stage2_mean']['auc']:.2f}",
    "AUC_S2_V3": f"{v3['both screened']['stage2_mean'][0]:.2f}",
    "SH_IN": f"{density(4, 20):.1f}", "SH_OUT": f"{density(-10, -4):.1f}",
    "Z_MED": f"{pn['z_median']:.2f}",
    "PZ_D_PCT": f"{100 * pn['n_ls_photo'] / (pn['n_ls_photo'] + pn['n_own_photo']):.0f}",
    "LOW": num["lowest"]["name"].split()[1], "HIGH": num["highest"]["name"].split()[1],
    "HALF": f"{ce['astrometry']['half_within']:.2f}",
    "REFUSED": f"{100 * ce['astrometry']['refused']:.0f}",
    "NEAR_PUR": f"{near['a: within 20']['purity']:.0f}", "FAR_PUR": f"{near['c: beyond 30']['purity']:.0f}",
    "S02_O": f"{above['original']['share_of_part_above']:.0f}", "S02_RE": f"{others:.0f}",
    "S02_ALL": f"{pp['score_share_all']['0.2']:.0f}",
    "S02_LO": f"{np.floor(lo02):.0f}", "S02_HI": f"{np.ceil(hi02):.0f}",
    "S02_LO1": f"{lo02:.1f}", "S02_HI1": f"{hi02:.1f}",
    "ZENODO": ZENODO,
    "TABLE_ROWS": "\n".join(rows),
}


def render(tpl, keys):
    used = set(re.findall(r"\{\{(\w+)\}\}", tpl))
    missing = used - set(keys)
    if missing:
        raise SystemExit(f"keys with no value: {sorted(missing)}")
    out = re.sub(r"\{\{(\w+)\}\}", lambda m: keys[m.group(1)], tpl)
    if "{{" in out:
        raise SystemExit("a placeholder survived the render")
    return out, sorted(set(keys) - used)


def unit(ra, dec):
    a, b = np.radians(ra), np.radians(dec)
    return np.column_stack([np.cos(b) * np.cos(a), np.cos(b) * np.sin(a), np.sin(b)])


def rules_hold():
    """The Note says rules (4) and (5) and the merge hold for every entry.
    Count the entries that break them, in the merged catalogue: another entry
    within max(6", r_eff), a DESI DR1 ZWARN = 0 spectrum within 2", three or
    more other entries within 1.5', or the excluded region. Planted: two
    entries 5.5" apart, entries moved 0.3" from DESI spectra, and four entries
    put within 1' of each other, must be caught."""
    from scipy.spatial import cKDTree
    sys.path.insert(0, "/home/ben/star/sweeps")
    import badregion
    d = np.load(f"{C}/desi_dr1_zwarn0.npz")
    tree = cKDTree(unit(d["ra"], d["dec"]))
    chord = 2 * np.sin(np.radians(2 / 3600) / 2)

    def desi(ra, dec):
        return np.isfinite(tree.query(unit(ra, dec), distance_upper_bound=chord)[0])

    if desi(d["ra"][:50] + 0.3 / 3600 / np.cos(np.radians(d["dec"][:50])), d["dec"][:50]).sum() != 50:
        raise SystemExit("planted DESI matches not caught")
    pr, pd = 10.0 + np.array([0, 0.01, 0.0, 0.01]), -40.0 + np.array([0, 0, 0.01, 0.01])
    if not badregion.flag(pr, pd).all():
        raise SystemExit("planted clump not caught")
    def merged(ra, dec, r):
        """Entries with another within max(6", the larger r_eff)."""
        x = unit(ra, dec)
        p = cKDTree(x).query_pairs(2 * np.sin(np.radians(max(6, r.max()) / 3600) / 2), output_type="ndarray")
        sep = np.degrees(2 * np.arcsin(np.linalg.norm(x[p[:, 0]] - x[p[:, 1]], axis=1) / 2)) * 3600
        p = p[sep < np.maximum(6, np.maximum(r[p[:, 0]], r[p[:, 1]]))]
        bad = np.zeros(len(ra), bool)
        bad[p.ravel()] = True
        return bad

    pr2 = np.array([10.0, 10.0 + 5.5 / 3600 / np.cos(np.radians(-40))])
    if merged(pr2, np.array([-40.0, -40.0]), np.array([3.0, 3.0])).sum() != 2:
        raise SystemExit("planted pair 5.5\" apart not caught")
    out = {}
    for name, bad in (("merge (another within max(6\", r_eff))", merged(z["ra"], z["dec"], z["shape_r"])),
                      ("DESI within 2\"", desi(z["ra"], z["dec"])),
                      ("clump (3+ others within 1.5')", badregion.flag(z["ra"], z["dec"]) & ~badregion.in_blacklist(z["ra"], z["dec"])),
                      ("excluded region", badregion.in_blacklist(z["ra"], z["dec"]))):
        out[name] = {s: int(bad[z["sub"] == s].sum()) for s in ("original", "recovered", "extension")}
    return out


def words(path):
    r = subprocess.run([TEXCOUNT, "-1", "-sum", "-incbib", "-merge", path],
                       capture_output=True, text=True)
    return int(r.stdout.strip().split()[0])


def main():
    tpl = open(TPL).read()
    # planted: a template key with no value must stop the render
    try:
        render(tpl + "{{NO_SUCH_KEY}}", KEYS)
        raise SystemExit("planted missing key was not caught")
    except SystemExit as e:
        if "NO_SUCH_KEY" not in str(e):
            raise
    out, unused = render(tpl, KEYS)
    if unused:
        print(f"keys not used by the template: {unused}")
    open(OUT, "w").write(out)
    if ce.get("draft"):
        print(f"NOTE: astrometry and scores are from a draft run: {ce['missing']:,} entries not yet "
              f"profiled, {ce['scored']:,} of {num['entries']:,} scored")
    if ZENODO.startswith(r"\textbf"):
        print("NOTE: Zenodo DOI not yet in the text")
    for name, by in rules_hold().items():
        if sum(by.values()):
            print(f"CLAIM FALSE: rule says none, but {sum(by.values())} entries break it: {name} {by}")
    ab = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", out, re.S).group(1)
    tmp = f"{HERE}/.abstract.tex"
    open(tmp, "w").write(ab)
    na, nt = words(tmp), words(OUT)
    os.remove(tmp)
    print(f"words: {nt} of 1,500 (texcount, references included); abstract {na} of 150")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
