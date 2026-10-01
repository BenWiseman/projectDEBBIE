"""Does a second-pass classifier on human labels have signal? Measure before building.

Joins the human verdicts (purity sample + shortlist) to the distilled source
table and cross-validates a gradient-boosted classifier on catalogue features
that were never part of the selection. The question is only whether the signal
exists; the production design is a separate matter.

The join is vectorised on an integer position key. An earlier version built a
33.5-million-entry Python dict, reached 11.7 GB and was killed.
"""
from __future__ import annotations

import csv
import io
import json
import sys

import numpy as np

SRC = "/home/ben/legacy/resolved_sources.npz"
PURITY = "/home/ben/Downloads/purity_verdicts_v2.jsonl"
SHORTLIST = "/home/ben/star/catalogue_debbie/shortlist_verdicts_2026-09-30.csv"

NAMES = ["shape_r", "shape_r_ivar", "r_signif", "sn_g", "sn_r", "dchisq_diff",
         "dchisq_ratio", "sersic", "e1", "e2", "ellip", "fracflux_g", "fracflux_r",
         "fracmasked_g", "fracin_g", "ebv", "nobs_g", "nobs_r", "g_r", "r_z", "r_w1"]


def poskey(ra, dec):
    """int64 key from a position at 5 decimal places."""
    return (np.round(np.asarray(ra, float) * 1e5).astype(np.int64) * 200_000_00
            + np.round((np.asarray(dec, float) + 90.0) * 1e5).astype(np.int64))


def load_labels():
    ra, dec, v, w, src = [], [], [], [], []
    raw = io.open(PURITY, encoding="utf-8").read()
    for line in raw.replace(chr(92) + "n", chr(10)).split(chr(10)):
        if line.strip():
            d = json.loads(line)
            ra.append(d["ra"]); dec.append(d["dec"]); v.append(d["v"])
            w.append(float(d.get("w", 1.0))); src.append("purity")
    for row in csv.DictReader(io.open(SHORTLIST, encoding="utf-8")):
        ra.append(float(row["ra_deg"])); dec.append(float(row["dec_deg"]))
        v.append(row["verdict"]); w.append(1.0); src.append("shortlist")
    return (np.array(ra), np.array(dec), np.array(v), np.array(w), np.array(src))


def features(r):
    g = np.asarray(r["flux_g"], float); rr = np.asarray(r["flux_r"], float)
    z = np.asarray(r["flux_z"], float); w1 = np.asarray(r["flux_w1"], float)
    ig = np.asarray(r["ivar_g"], float); ir = np.asarray(r["ivar_r"], float)
    sr = np.asarray(r["shape_r"], float); sri = np.asarray(r["shape_r_ivar"], float)
    dp = np.asarray(r["dchisq_psf"], float); ds = np.asarray(r["dchisq_ser"], float)
    e1 = np.asarray(r["e1"], float); e2 = np.asarray(r["e2"], float)
    with np.errstate(all="ignore"):
        X = np.column_stack([
            sr, sri, sr * np.sqrt(np.maximum(sri, 0)),
            g * np.sqrt(np.maximum(ig, 0)), rr * np.sqrt(np.maximum(ir, 0)),
            ds - dp, np.where(dp != 0, ds / dp, 0.0),
            np.asarray(r["sersic"], float), e1, e2, np.hypot(e1, e2),
            np.asarray(r["fracflux_g"], float), np.asarray(r["fracflux_r"], float),
            np.asarray(r["fracmasked_g"], float), np.asarray(r["fracin_g"], float),
            np.asarray(r["ebv"], float),
            np.asarray(r["nobs_g"], float), np.asarray(r["nobs_r"], float),
            -2.5 * np.log10(np.maximum(g, 1e-3)) + 2.5 * np.log10(np.maximum(rr, 1e-3)),
            -2.5 * np.log10(np.maximum(rr, 1e-3)) + 2.5 * np.log10(np.maximum(z, 1e-3)),
            -2.5 * np.log10(np.maximum(rr, 1e-3)) + 2.5 * np.log10(np.maximum(w1, 1e-3)),
        ])
    return np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)


def main():
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_auc_score
    import xgboost as xgb

    R = np.load(SRC)["rows"]
    print(f"distilled sources {len(R):,}")
    cat_key = poskey(R["ra"], R["dec"])
    order = np.argsort(cat_key, kind="stable")
    sorted_key = cat_key[order]

    lra, ldec, lv, lw, lsrc = load_labels()
    lk = poskey(lra, ldec)
    pos = np.searchsorted(sorted_key, lk)
    pos = np.clip(pos, 0, len(sorted_key) - 1)
    ok = sorted_key[pos] == lk
    idx = order[pos[ok]]
    lv, lw, lsrc = lv[ok], lw[ok], lsrc[ok]
    print(f"labels {len(lk):,} -> joined {ok.sum():,} ({100 * ok.mean():.1f}%)")

    X_all = features(R[idx])

    for pop in ("purity", "both"):
        m = np.ones(len(lv), bool) if pop == "both" else (lsrc == "purity")
        d = m & np.isin(lv, ["yes", "no"])
        X, y, w = X_all[d], (lv[d] == "yes").astype(int), lw[d]
        oof = np.zeros(len(y))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=7).split(X, y):
            clf = xgb.XGBClassifier(n_estimators=400, max_depth=4, learning_rate=0.05,
                                    subsample=0.8, colsample_bytree=0.8,
                                    eval_metric="logloss", verbosity=0)
            clf.fit(X[tr], y[tr])
            oof[te] = clf.predict_proba(X[te])[:, 1]
        print(f"\n{pop:9s} n={len(y):4d}  base rate {y.mean():.3f}   "
              f"AUC {roc_auc_score(y, oof):.4f}   "
              f"weighted AUC {roc_auc_score(y, oof, sample_weight=w):.4f}")
        if pop == "purity":
            # what a threshold would actually do to catalogue purity, weighted
            print(f"   {'keep top':>9s} {'kept':>6s} {'weighted purity':>16s}")
            for frac in (1.0, 0.9, 0.8, 0.7, 0.6, 0.5):
                t = np.quantile(oof, 1 - frac)
                k = oof >= t
                p = w[k & (y == 1)].sum() / w[k].sum()
                print(f"   {100*frac:8.0f}% {k.sum():6d} {100*p:15.2f}%")
            np.savez("/home/ben/star/output/secondpass_probe.npz",
                     oof=oof, y=y, w=w)
    return 0


if __name__ == "__main__":
    sys.exit(main())
