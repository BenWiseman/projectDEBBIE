"""Measure the selection function empirically, against real catalogued galaxies.

A catalogue paper needs to say what fraction of the galaxies that exist in the
footprint actually made it in, as a function of size and surface brightness.
Without that, 243,781 is a number with no denominator.

Method
------
Injection-recovery into images is impossible now: the sweeps are distilled and
the survey imaging service is unreachable. But there is a better empirical
measurement available, and it uses real objects rather than synthetic ones.

Take every galaxy that a published catalogue already knows about inside our
footprint. Ask, for each: did the DR10 sweeps detect it, and would our selection
have kept it? That is end-to-end recovery of real galaxies with real profiles,
real neighbours and real imaging defects -- everything injection has to
approximate.

Two numbers come out, and they must not be conflated:

  DETECTION   fraction of known galaxies present in the sweep catalogue at all.
              This is the survey's completeness, not ours.
  SELECTION   fraction of those the cuts keep. This is ours, and it is the one
              a reader needs to interpret the catalogue.

The honest total is the product, reported per (r_eff, mu) bin.

What this cannot measure
------------------------
Galaxies nobody has catalogued anywhere. If a whole population is missing from
every published catalogue AND from the sweeps, this method is blind to it by
construction. That is stated rather than hidden: the measurement is completeness
relative to known galaxies, which is a lower bound on the true selection
function only if the known population is representative.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, "/home/ben/star/src")
sys.path.insert(0, "/home/ben/legacy")

LEGACY = "/home/ben/legacy"
OUT = os.path.join(LEGACY, "completeness.json")
EXTRACT = os.path.join(LEGACY, "resolved_sources.npz")

MIN_R, MIN_MU, MAX_MU = 3.0, 24.0, 26.5
RESOLVED = {"EXP", "DEV", "SER", "REX"}
MATCH_ARCSEC = 3.0


def mag(f):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(f > 0, 22.5 - 2.5 * np.log10(np.maximum(f, 1e-12)), np.nan)


def main():
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from scipy.spatial import cKDTree
    from starfinder.localcat import _unit

    print("loading the extract ...", flush=True)
    rows = np.load(EXTRACT, allow_pickle=False)["rows"]
    print(f"  {len(rows):,} resolved sources", flush=True)

    # known galaxies inside the footprint, from the same index used to screen
    print("loading known galaxies ...", flush=True)
    z = np.load(os.path.join(LEGACY, "catalogues", "named_index.npz"),
                allow_pickle=True)
    kra, kdec = z["ra"], z["dec"]

    # restrict to our footprint: |b| >= 30 and the RA/Dec span the sweeps cover
    b = np.abs(SkyCoord(kra * u.deg, kdec * u.deg).galactic.b.deg)
    inside = b >= 30.0
    tree_src = cKDTree(_unit(rows["ra"], rows["dec"]))
    print(f"  {inside.sum():,} known objects at |b|>=30", flush=True)

    rng = np.random.default_rng(0)
    idx = np.flatnonzero(inside)
    if len(idx) > 400_000:
        idx = rng.choice(idx, 400_000, replace=False)
    kra, kdec = kra[idx], kdec[idx]

    # RESTRICT TO THE OBSERVED FOOTPRINT.
    #
    # The index is all-sky; the sweeps cover only DR10 south. Comparing the two
    # directly asks "what fraction of all known galaxies lie in DR10 south",
    # which is a footprint ratio, not a completeness. It gave 8.4% and would
    # have been a badly wrong number in a paper.
    #
    # A position is inside the footprint if the sweeps contain ANY source
    # within 1 arcmin -- sources are dense enough (33.5M over the footprint)
    # that this marks observed sky, independently of whether this particular
    # galaxy was detected.
    foot_chord = 2.0 * np.sin(np.radians(60.0 / 3600.0) / 2.0)
    dfoot, _ = tree_src.query(_unit(kra, kdec), k=1,
                              distance_upper_bound=foot_chord)
    in_footprint = np.isfinite(dfoot) & (dfoot <= foot_chord)
    print(f"  {in_footprint.sum():,} of those fall inside the observed "
          f"footprint ({in_footprint.mean()*100:.1f}%)", flush=True)
    kra, kdec = kra[in_footprint], kdec[in_footprint]
    if len(kra) < 1000:
        print("too few known galaxies inside the footprint to measure")
        return

    # 1. DETECTION: is the known galaxy in the sweep extract at all?
    chord = 2.0 * np.sin(np.radians(MATCH_ARCSEC / 3600.0) / 2.0)
    d, j = tree_src.query(_unit(kra, kdec), k=1, distance_upper_bound=chord)
    detected = np.isfinite(d) & (d <= chord)
    j_ok = np.where(detected, j, 0)
    # This is NOT detection completeness, and must not be labelled as such.
    # The extract keeps only resolved sources with r_eff >= 2", and most
    # catalogued galaxies are smaller than that, so this number conflates "DR10
    # detected it" with "it is larger than 2 arcsec". The two cannot be
    # separated now that the raw sweeps are distilled. Reported for what it is.
    print(f"\nIN EXTRACT: {detected.sum():,}/{len(kra):,} "
          f"({detected.mean()*100:.1f}%) of known galaxies inside the footprint "
          f"are present as resolved sources with r_eff >= 2 arcsec.",
          flush=True)
    print("  NOTE: this is not detection completeness -- it folds in the 2 "
          "arcsec size cut of the extract, and most catalogued galaxies are "
          "smaller than that.", flush=True)

    # 2. SELECTION: of those, how many pass our cuts?
    sub = rows[j_ok]
    g = mag(sub["flux_g"])
    re = sub["shape_r"].astype(np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        mu = g + 2.5 * np.log10(2 * np.pi * np.maximum(re, 1e-3) ** 2)
    passes = (np.isin(sub["type"], list(RESOLVED)) & (re >= MIN_R)
              & (sub["nobs_g"] >= 2) & (sub["nobs_r"] >= 2)
              & (sub["maskbits"] == 0)
              & (sub["flux_g"] > 0) & (sub["flux_r"] > 0) & (sub["flux_z"] > 0)
              & (mu >= MIN_MU) & (mu <= MAX_MU) & np.isfinite(mu))
    passes &= detected
    print(f"SELECTION: {passes.sum():,}/{detected.sum():,} "
          f"({passes.sum()/max(detected.sum(),1)*100:.1f}%) of detected "
          f"galaxies pass the cuts", flush=True)
    print(f"\nThe number a reader needs is the SELECTION FUNCTION below: given "
          f"a galaxy present in the sweeps at r_eff >= 3 arcsec, what fraction "
          f"of the cuts keep it.", flush=True)

    # 3. the selection function, binned -- what a reader actually needs
    print("\nselection fraction by r_eff and mu (detected galaxies only):",
          flush=True)
    grid = {}
    r_edges = [3, 4, 5, 6, 8, 10, 15, 100]
    m_edges = [22, 23, 24, 25, 26, 26.5, 30]
    print(f"  {'r_eff':>10s} " + " ".join(f"{a:>5.1f}" for a in m_edges[:-1]))
    for i in range(len(r_edges) - 1):
        line = f"  {r_edges[i]:4.0f}-{r_edges[i+1]:<5.0f}"
        for k in range(len(m_edges) - 1):
            m = (detected & (re >= r_edges[i]) & (re < r_edges[i + 1])
                 & (mu >= m_edges[k]) & (mu < m_edges[k + 1]))
            if m.sum() < 30:
                line += "     ."
                continue
            frac = passes[m].mean()
            grid[f"{r_edges[i]}-{r_edges[i+1]}_{m_edges[k]}-{m_edges[k+1]}"] = \
                {"n": int(m.sum()), "frac": float(frac)}
            line += f" {frac*100:5.0f}"
        print(line, flush=True)

    json.dump({"n_known_sampled": int(len(kra)),
               "in_extract_frac": float(detected.mean()),
               "in_extract_note": "NOT detection completeness: folds in the "
                                  "extract's r_eff >= 2 arcsec cut",
               "selection_frac": float(passes.sum() / max(detected.sum(), 1)),

               "grid": grid,
               "cuts": {"min_r": MIN_R, "min_mu": MIN_MU, "max_mu": MAX_MU},
               "note": "completeness relative to KNOWN galaxies; blind to any "
                       "population absent from both published catalogues and "
                       "the sweeps"},
              open(OUT, "w"), indent=2)
    print(f"\nsaved {OUT}")


if __name__ == "__main__":
    main()
