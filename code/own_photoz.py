"""Estimate our own photometric redshifts, where no catalogue covers the sky.

Public photo-z gets us 55.8% and stops. The gap is not depth -- matched and
unmatched objects have identical magnitudes, sizes and surface brightness, and
the match rate is ~55% at every magnitude. It is footprint: VII/292 covers
85.9% of our objects at RA 0-30 and 6.6% at RA 270-300.

So no amount of further querying helps. But we measured g, r, z, W1 and W2 for
100% of the catalogue, and 135,140 objects have a redshift from a public
catalogue. That is a training set for a straightforward broadband photo-z.

What this is and is not
-----------------------
It is a colour-based estimate from five bands, trained on objects that a public
photo-z already covered. It will be worse than a dedicated photo-z pipeline
with more bands and a spectroscopic training set, and it inherits any bias in
its training labels -- which are themselves photometric for all but ~2,000
objects. Errors compound.

It is reported as its own source ("own-photoz") and never merged silently with
catalogue redshifts. The honest use is population statistics, not per-object
distances.

Validation is on a HELD-OUT SKY REGION, not random rows: photo-z errors are
spatially correlated through calibration and extinction, so a random split
flatters the result.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, "/home/ben/legacy")

LEGACY = "/home/ben/legacy"
OUT = os.path.join(LEGACY, "own_photoz.npz")


def build():
    cat = np.load(os.path.join(LEGACY, "VERIFIED_NOVEL.npy"))
    m = np.load(os.path.join(LEGACY, "measurements.npz"), allow_pickle=True)
    r = np.load(os.path.join(LEGACY, "redshifts_all.npz"), allow_pickle=True)

    # Pull every array out of the npz ONCE. Indexing `m["g"][i]` inside a loop
    # re-decompresses the whole column on every access -- six columns across
    # 243,727 objects made this take longer than a 28-minute timeout.
    mg, mr, mz = m["g"], m["r"], m["z"]
    mw1 = m["w1"] if "w1" in m.files else np.full(len(mg), np.nan)
    m_ra, m_dec = m["ra"], m["dec"]
    r_ra, r_dec, r_z = r["ra"], r["dec"], r["z"]

    mk = {f"{a:.6f}_{b:.6f}": i for i, (a, b) in enumerate(zip(m_ra, m_dec))}
    zk = {f"{a:.6f}_{b:.6f}": zz for a, b, zz in zip(r_ra, r_dec, r_z)}

    print(f"  indexing {len(cat):,} objects ...", flush=True)
    idx, keep = [], []
    for j, o in enumerate(cat):
        i = mk.get(f"{o['ra']:.6f}_{o['dec']:.6f}")
        if i is not None:
            idx.append(i)
            keep.append(j)
    idx = np.array(idx)
    keep = np.array(keep)
    print(f"  {len(idx):,} matched to measurements", flush=True)

    # vectorised feature build
    X = np.column_stack([mg[idx], mr[idx], mz[idx],
                         mg[idx] - mr[idx], mr[idx] - mz[idx], mw1[idx],
                         cat["shape_r"][keep], cat["mu_g"][keep]]).astype(np.float32)
    z = np.array([zk.get(f"{a:.6f}_{b:.6f}", np.nan)
                  for a, b in zip(cat["ra"][keep], cat["dec"][keep])])
    return X, z, cat["ra"][keep], cat["dec"][keep]


def main():
    X, z, ra, dec = build()
    known = np.isfinite(z) & (z > 0) & np.isfinite(X[:, :5]).all(axis=1)
    print(f"{len(X):,} objects, {known.sum():,} with a training redshift")
    if known.sum() < 5000:
        print("too few labels")
        return

    # held-out SKY, not random rows: photo-z errors are spatially correlated
    box = (np.floor(ra / 15).astype(int) * 100 + np.floor(dec / 15).astype(int))
    ub = np.unique(box[known])
    rng = np.random.default_rng(0)
    rng.shuffle(ub)
    test_boxes = set(ub[: max(1, len(ub) // 4)].tolist())
    is_test = np.array([b in test_boxes for b in box])
    tr = known & ~is_test
    te = known & is_test
    print(f"  train {tr.sum():,} / test {te.sum():,} over {len(ub)} sky boxes")

    from sklearn.ensemble import HistGradientBoostingRegressor
    mdl = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.08,
                                        early_stopping=True, random_state=0)
    mdl.fit(np.nan_to_num(X[tr], nan=-99), z[tr])
    pred = mdl.predict(np.nan_to_num(X[te], nan=-99))

    dz = (pred - z[te]) / (1 + z[te])
    sigma = 1.4826 * np.median(np.abs(dz - np.median(dz)))
    outlier = np.mean(np.abs(dz) > 0.15)
    print(f"\nheld-out sky:")
    print(f"  sigma_NMAD    {sigma:.4f}")
    print(f"  outlier rate  {outlier*100:.1f}%  (|dz|/(1+z) > 0.15)")
    print(f"  bias          {np.median(dz):+.4f}")

    # a control: shuffle the labels, the metric must collapse
    ctl = HistGradientBoostingRegressor(max_iter=100, early_stopping=True,
                                        random_state=0)
    ctl.fit(np.nan_to_num(X[tr], nan=-99), rng.permutation(z[tr]))
    cdz = (ctl.predict(np.nan_to_num(X[te], nan=-99)) - z[te]) / (1 + z[te])
    csig = 1.4826 * np.median(np.abs(cdz - np.median(cdz)))
    print(f"  shuffled control sigma_NMAD {csig:.4f}  (must be much worse)")

    if sigma > 0.15 or sigma > csig * 0.8:
        print("\nNOT USABLE: the model is not meaningfully better than chance.")
        return

    # predict for everything without a catalogue redshift
    mdl.fit(np.nan_to_num(X[known], nan=-99), z[known])
    need = ~known & np.isfinite(X[:, :5]).all(axis=1)
    zown = np.full(len(X), np.nan)
    zown[need] = mdl.predict(np.nan_to_num(X[need], nan=-99))
    zown[zown < 0] = np.nan
    np.savez_compressed(OUT, ra=ra, dec=dec, z_own=zown,
                        sigma_nmad=sigma, outlier_rate=outlier,
                        n_train=int(known.sum()))
    print(f"\npredicted for {np.isfinite(zown).sum():,} objects without one")
    print(f"  combined coverage would be "
          f"{(known | np.isfinite(zown)).sum():,}/{len(X):,} "
          f"({(known | np.isfinite(zown)).mean()*100:.1f}%)")
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
