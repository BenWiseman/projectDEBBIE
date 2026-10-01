"""Catalogue purity measured separately in each (r_eff, mu) regime.

The headline 70.3% averages over a catalogue that is not homogeneous. Split on
surface brightness and the two halves are 95.6% and 50.0%, which is a different
paper: the bright half is a clean catalogue nine times the size of SMUDGes,
and the faint half is a flagged candidate list.

The cut is on mu_eff, a released column, not on a model score. That matters:
the selection function of Section 5.5 is measured on exactly these axes and
stays analytic, there is no train/test split to worry about, and nothing in the
headline product depends on a classifier.

    .venv/bin/python sweeps/purity_by_regime.py [verdicts.jsonl]
"""
from __future__ import annotations

import io
import json
import math
import sys

import numpy as np

CAT = "/home/ben/star/catalogue_debbie/v1.2_main.npy"
DEFAULT_VERDICTS = "/home/ben/Downloads/purity_verdicts_v2.jsonl"
Z = 1.959964
MU_SPLIT = 25.5
R_SPLIT = 6.0


def wilson(p, n):
    if n <= 0:
        return (0.0, 0.0)
    z2 = Z * Z
    den = 1 + z2 / n
    c = (p + z2 / (2 * n)) / den
    h = (Z / den) * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n))
    return (max(0.0, c - h), min(1.0, c + h))


def weighted(rows):
    """Weighted purity and the Wilson interval at the effective sample size."""
    wy = sum(r["w"] for r in rows if r["v"] == "yes")
    wn = sum(r["w"] for r in rows if r["v"] == "no")
    sw, sw2 = wy + wn, sum(r["w"] ** 2 for r in rows)
    if sw <= 0:
        return None
    p = wy / sw
    ne = sw * sw / sw2 if sw2 else 0.0
    return p, wilson(p, ne), ne, sw


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_VERDICTS
    raw = io.open(path, encoding="utf-8").read()
    rows = [json.loads(l) for l in raw.replace(chr(92) + "n", chr(10)).split(chr(10)) if l.strip()]
    dec = [r for r in rows if r["v"] in ("yes", "no")]
    und = len(rows) - len(dec)
    print(f"verdicts {len(rows):,}  decided {len(dec):,}  undecided {und:,}\n")

    cat = np.load(CAT, allow_pickle=True)
    cr = np.asarray(cat["shape_r"], float)
    cmu = np.asarray(cat["mu_g"], float)
    n_cat = len(cr)

    regimes = [
        (f"r>={R_SPLIT:.0f}, mu<{MU_SPLIT}", lambda x: x["r"] >= R_SPLIT and x["mu"] < MU_SPLIT,
         (cr >= R_SPLIT) & (cmu < MU_SPLIT)),
        (f"r<{R_SPLIT:.0f},  mu<{MU_SPLIT}", lambda x: x["r"] < R_SPLIT and x["mu"] < MU_SPLIT,
         (cr < R_SPLIT) & (cmu < MU_SPLIT)),
        (f"r>={R_SPLIT:.0f}, mu>={MU_SPLIT}", lambda x: x["r"] >= R_SPLIT and x["mu"] >= MU_SPLIT,
         (cr >= R_SPLIT) & (cmu >= MU_SPLIT)),
        (f"r<{R_SPLIT:.0f},  mu>={MU_SPLIT}", lambda x: x["r"] < R_SPLIT and x["mu"] >= MU_SPLIT,
         (cr < R_SPLIT) & (cmu >= MU_SPLIT)),
    ]

    print(f"{'regime':<20s} {'judged':>7s} {'in catalogue':>13s} {'purity':>8s} {'95% CI':>16s}")
    out = {}
    for label, pred, cmask in regimes:
        sub = [r for r in dec if pred(r)]
        res = weighted(sub)
        if not res:
            continue
        p, (lo, hi), ne, _ = res
        print(f"{label:<20s} {len(sub):7d} {int(cmask.sum()):13,d} "
              f"{100 * p:7.1f}% [{100 * lo:5.1f},{100 * hi:5.1f}]")
        out[label] = {"judged": len(sub), "n_catalogue": int(cmask.sum()),
                      "purity": round(p, 4), "ci": [round(lo, 4), round(hi, 4)],
                      "n_eff": round(ne, 1)}

    allres = weighted(dec)
    print(f"{'ALL':<20s} {len(dec):7d} {n_cat:13,d} {100 * allres[0]:7.1f}% "
          f"[{100 * allres[1][0]:5.1f},{100 * allres[1][1]:5.1f}]")

    # the headline split, on surface brightness alone
    print()
    for label, pred, cmask in (("bright half", lambda x: x["mu"] < MU_SPLIT, cmu < MU_SPLIT),
                               ("faint half", lambda x: x["mu"] >= MU_SPLIT, cmu >= MU_SPLIT)):
        sub = [r for r in dec if pred(r)]
        p, (lo, hi), ne, _ = weighted(sub)
        n = int(cmask.sum())
        print(f"{label:<12s} mu {'<' if 'bright' in label else '>='}{MU_SPLIT}: "
              f"{n:,} entries ({100 * n / n_cat:.1f}% of the catalogue), "
              f"purity {100 * p:.1f}% [{100 * lo:.1f}, {100 * hi:.1f}]")
        out[label] = {"judged": len(sub), "n_catalogue": n, "purity": round(p, 4),
                      "ci": [round(lo, 4), round(hi, 4)], "n_eff": round(ne, 1)}

    out["ALL"] = {"judged": len(dec), "n_catalogue": n_cat,
                  "purity": round(allres[0], 4),
                  "ci": [round(allres[1][0], 4), round(allres[1][1], 4)],
                  "n_eff": round(allres[2], 1), "undecided": und}
    dest = "/home/ben/star/catalogue_debbie/purity_by_regime.json"
    io.open(dest, "w", encoding="utf-8").write(json.dumps(out, indent=2) + "\n")
    print(f"\nwrote {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
