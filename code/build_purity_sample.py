"""Draw the purity sample, stratified, from v1.2.

Two defects in the first draw, both found by Ben:

  1. It was drawn from v1.1, before DESI DR1 was screened. 9.15% of it was
     objects DESI has already identified spectroscopically -- confirmed
     galaxies that are no longer in the catalogue. Leaving them in biases the
     purity estimate upward, which is the direction that flatters us.

  2. It was a uniform draw. Uniform is the correct weighting for catalogue
     purity -- every entry counts once -- but it spends 85% of the effort at
     r_eff 3-4 arcsec and lands about ONE object in each of the largest cells,
     so it says nothing about the size range worth caring about.

The fix is stratification with weights, not abandoning uniformity: sample the
(r_eff, mu) plane in cells, give every object the weight of the catalogue
entries it stands for, and reconstruct the unbiased total by weighting.

Allocation is hybrid, 70% proportional and 30% equal, chosen by measuring the
trade rather than by taste:

    scheme                 n_eff   smallest cell
    equal per cell           841              71
    proportional (uniform)  1998               1
    sqrt-proportional       1655              12
    hybrid 70% proportional 1776              22   <-- chosen

An 11% loss on the headline number buys a factor of 22 in the cells that were
otherwise empty. The error on the weighted purity uses the effective sample
size, not a naive binomial over the object count, because the weights are
unequal.
"""
from __future__ import annotations

import json
import os

import numpy as np

CAT = "/home/ben/star/catalogue_debbie/v1.2_main.npy"
OUT = "/home/ben/legacy/inspect_sample_v2.json"
R_EDGES = [3, 3.5, 4, 5, 6, 8, 10, 15.1]
MU_EDGES = [24, 25, 25.5, 26, 26.51]
TARGET = 2000
PROP_FRACTION = 0.7
SEED = 20260905


def main():
    m = np.load(CAT)
    cells = []
    for i in range(len(R_EDGES) - 1):
        for j in range(len(MU_EDGES) - 1):
            sel = np.flatnonzero(
                (m["shape_r"] >= R_EDGES[i]) & (m["shape_r"] < R_EDGES[i + 1])
                & (m["mu_g"] >= MU_EDGES[j]) & (m["mu_g"] < MU_EDGES[j + 1]))
            if len(sel):
                cells.append((f"r{R_EDGES[i]}-{R_EDGES[i+1]} "
                              f"mu{MU_EDGES[j]}-{MU_EDGES[j+1]}", sel))

    n_cell = np.array([len(s) for _, s in cells], dtype=float)
    alloc = (PROP_FRACTION * TARGET * n_cell / n_cell.sum()
             + (1 - PROP_FRACTION) * TARGET / len(cells))
    alloc = np.minimum(np.maximum(np.round(alloc), 1), n_cell).astype(int)

    rng = np.random.default_rng(SEED)
    objs = []
    for (name, sel), n in zip(cells, alloc):
        pick = rng.choice(sel, n, replace=False)
        w = len(sel) / n                      # entries this object stands for
        for i in pick:
            o = m[i]
            objs.append({
                "ra": float(o["ra"]), "dec": float(o["dec"]),
                "r": float(o["shape_r"]), "mu": float(o["mu_g"]),
                "w": float(w), "cell": name,
                "img": (f"candidate_images/ra{int(o['ra']):03d}/"
                        f"{o['ra']:.5f}{o['dec']:+.5f}_75.jpg"),
            })
    rng.shuffle(objs)

    total_w = sum(o["w"] for o in objs)
    n_eff = total_w ** 2 / sum(o["w"] ** 2 for o in objs)
    # the weights must reconstruct the catalogue, or the estimate is not unbiased
    assert abs(total_w - len(m)) < 0.02 * len(m), (total_w, len(m))
    missing = [o for o in objs if not os.path.exists("/home/ben/legacy/" + o["img"])]
    print(f"{len(objs):,} objects over {len(cells)} cells")
    print(f"  weights sum {total_w:,.0f} vs catalogue {len(m):,}")
    print(f"  effective n for the weighted headline: {n_eff:,.0f}")
    print(f"  missing cutouts: {len(missing)}")

    json.dump({"seed": SEED,
               "scheme": f"hybrid {int(PROP_FRACTION*100)}% proportional",
               "source": os.path.basename(CAT), "n": len(objs),
               "n_eff": float(n_eff), "objects": objs}, open(OUT, "w"))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
