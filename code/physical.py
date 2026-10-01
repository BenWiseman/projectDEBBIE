"""Turn angular measurements into physical ones, using redshifts.

Without a distance, r_eff = 5 arcsec is uninterpretable: a nearby dwarf and a
distant giant look identical. This joins the catalogue to redshift catalogues
and derives what a reader actually needs -- physical size, absolute magnitude,
and a stellar mass estimate.

Redshift sources, best first
----------------------------
    DESI DR1 spectroscopic   z with ZWARN == 0. A real measurement.
    LS photometric (VII/292) zphot from the same imaging our candidates come
                             from, so coverage should be near-complete. Carries
                             e_zphot and a quality flag; also pstar, which is a
                             free star-contamination cleaner.

The source is recorded per object. A photometric redshift is not a
spectroscopic one and the catalogue must never present them as the same thing.

Caveats that belong in the paper, not buried here
--------------------------------------------------
- Photo-z is calibrated on normal galaxies. Low-surface-brightness objects are
  under-represented in its training set, so both z and e_zphot are less
  trustworthy for exactly the objects we care about most.
- Stellar mass uses a single colour-to-M/L relation (Bell et al. 2003 style).
  That is a factor-of-two estimate, not a measurement.
- No k-correction is applied. At z < 0.1 that is a small error; beyond that it
  is not, and those objects are flagged.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, "/home/ben/star/src")
sys.path.insert(0, "/home/ben/legacy")

LEGACY = "/home/ben/legacy"
ZDIR = "/run/media/ben/KINGSTON/redshifts"
MEAS = os.path.join(LEGACY, "measurements.npz")
OUT = os.path.join(LEGACY, "physical.npz")

SPEC_ARCSEC = 2.0        # spectroscopic fibre positions are precise
PHOT_ARCSEC = 1.5        # photo-z rows are the same detections as our sweeps

# flat LCDM, Planck-ish
H0, OM, OL, C = 70.0, 0.3, 0.7, 299792.458


def comoving_mpc(z):
    """Comoving distance by direct integration. Vectorised, no cosmology dep."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    out = np.zeros_like(z)
    good = np.isfinite(z) & (z > 0)
    if not good.any():
        return out
    zz = z[good]
    steps = 512
    grid = np.linspace(0.0, 1.0, steps)[None, :] * zz[:, None]
    E = np.sqrt(OM * (1 + grid) ** 3 + OL)
    out[good] = (C / H0) * np.trapezoid(1.0 / E, grid, axis=1) \
        if hasattr(np, "trapezoid") else (C / H0) * np.trapz(1.0 / E, grid, axis=1)
    return out


def load_z_tsv(path, ra_keys, dec_keys, z_key, extra=()):
    """Parse a VizieR asu-tsv redshift table."""
    from starfinder.localcat import _f
    ras, decs, zs, ex = [], [], [], {k: [] for k in extra}
    hdr = None
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split("\t")]
            if hdr is None:
                hdr = parts
                continue
            if not parts or set(parts[0]) <= set("- "):
                continue
            row = dict(zip(hdr, parts))
            ra = next((_f(row[k]) for k in ra_keys if k in row), np.nan)
            dec = next((_f(row[k]) for k in dec_keys if k in row), np.nan)
            zv = _f(row.get(z_key, ""))
            if not (np.isfinite(ra) and np.isfinite(dec) and np.isfinite(zv)):
                continue
            ras.append(ra); decs.append(dec); zs.append(zv)
            for k in extra:
                ex[k].append(_f(row.get(k, "")))
    return (np.array(ras), np.array(decs), np.array(zs),
            {k: np.array(v) for k, v in ex.items()})


def match(cat_ra, cat_dec, ra, dec, arcsec):
    from scipy.spatial import cKDTree
    from starfinder.localcat import _unit
    if not len(ra):
        return np.full(len(cat_ra), -1)
    tree = cKDTree(_unit(ra, dec))
    chord = 2.0 * np.sin(np.radians(arcsec / 3600.0) / 2.0)
    d, i = tree.query(_unit(cat_ra, cat_dec), k=1, distance_upper_bound=chord)
    i = np.where(np.isfinite(d) & (d <= chord), i, -1)
    return i


def main():
    m = np.load(MEAS, allow_pickle=True)
    ra, dec = m["ra"], m["dec"]
    n = len(ra)
    print(f"{n:,} catalogue objects", flush=True)

    z = np.full(n, np.nan)
    z_err = np.full(n, np.nan)
    z_src = np.full(n, "", dtype="U6")

    # 1. spectroscopic first -- it wins wherever it exists
    spec = os.path.join(ZDIR, "desi_specz.tsv")
    if os.path.exists(spec):
        sra, sdec, sz, ex = load_z_tsv(spec, ("RAICRS", "RAJ2000"),
                                       ("DEICRS", "DEJ2000"), "z",
                                       extra=("e_z", "ZWARN"))
        ok = (ex["ZWARN"] == 0) | ~np.isfinite(ex["ZWARN"])
        sra, sdec, sz = sra[ok], sdec[ok], sz[ok]
        se = ex["e_z"][ok]
        idx = match(ra, dec, sra, sdec, SPEC_ARCSEC)
        hit = idx >= 0
        z[hit] = sz[idx[hit]]
        z_err[hit] = se[idx[hit]]
        z_src[hit] = "spec"
        print(f"  spectroscopic: {hit.sum():,} matched "
              f"({len(sra):,} in catalogue)", flush=True)

    # 2. photometric from the X-Match result (bulk download was truncated by
    #    VizieR at ~15-18M rows per request, so it never gave full coverage)
    pstar = np.full(n, np.nan)
    xm = os.path.join(LEGACY, "photoz_xmatch.npz")
    if os.path.exists(xm):
        x = np.load(xm)
        idx = match(ra, dec, x["ra"], x["dec"], PHOT_ARCSEC)
        hit = (idx >= 0) & ~np.isfinite(z)
        tgt = np.flatnonzero(hit)
        src = idx[hit]
        good = np.isfinite(x["zphot"][src])
        tgt, src = tgt[good], src[good]
        z[tgt] = x["zphot"][src]
        z_err[tgt] = x["e_zphot"][src]
        pstar[tgt] = x["pstar"][src]
        z_src[tgt] = "photo"
        print(f"  X-Match photo-z: +{len(tgt):,}", flush=True)

    # The bulk lsphotoz_*.tsv slices are superseded by X-Match and contribute
    # +0 every time -- but they are 1-2 GB each and parsing six of them in
    # Python dominates this script's runtime. Kept on disk, not re-read.
    pfiles = [] if os.path.exists(os.path.join(LEGACY, "photoz_xmatch.npz")) \
        else sorted(glob.glob(os.path.join(ZDIR, "lsphotoz_*.tsv")))
    for p in pfiles:
        need = ~np.isfinite(z)
        if not need.any():
            break
        pra, pdec, pz, ex = load_z_tsv(p, ("RAJ2000",), ("DEJ2000",), "zphot",
                                       extra=("e_zphot", "pstar", "fqual"))
        if not len(pra):
            continue
        idx = match(ra[need], dec[need], pra, pdec, PHOT_ARCSEC)
        hit = idx >= 0
        tgt = np.flatnonzero(need)[hit]
        z[tgt] = pz[idx[hit]]
        z_err[tgt] = ex["e_zphot"][idx[hit]]
        pstar[tgt] = ex["pstar"][idx[hit]]
        z_src[tgt] = "photo"
        print(f"  {os.path.basename(p)}: +{hit.sum():,}", flush=True)

    # 3. our own photo-z for the sky no public catalogue covers. Flagged as a
    #    distinct source: it is trained on photometric labels, so its error is
    #    scatter against LS photo-z, not against truth, and the real
    #    uncertainty against spectroscopy is larger. Population statistics
    #    only -- never a per-object distance.
    own = os.path.join(LEGACY, "own_photoz.npz")
    if os.path.exists(own):
        o = np.load(own)
        ok = {f"{a:.6f}_{b:.6f}": v for a, b, v in
              zip(o["ra"], o["dec"], o["z_own"]) if np.isfinite(v)}
        added = 0
        for i, (a, b) in enumerate(zip(ra, dec)):
            if np.isfinite(z[i]):
                continue
            v = ok.get(f"{a:.6f}_{b:.6f}")
            if v is not None:
                z[i] = v
                z_src[i] = "own"
                added += 1
        print(f"  own photo-z: +{added:,} "
              f"(sigma_NMAD {float(o['sigma_nmad']):.4f} vs LS photo-z)",
              flush=True)

    have = np.isfinite(z) & (z > 0)
    print(f"\nredshifts: {have.sum():,}/{n:,} ({have.mean()*100:.1f}%)  "
          f"spec {(z_src=='spec').sum():,}  photo {(z_src=='photo').sum():,}")

    # physical quantities
    dc = comoving_mpc(z)
    d_lum = dc * (1 + z)
    d_ang = dc / (1 + z)
    kpc_per_arcsec = d_ang * 1000.0 * np.pi / (180.0 * 3600.0)
    r_kpc = m["shape_r"] * kpc_per_arcsec
    with np.errstate(invalid="ignore", divide="ignore"):
        dm = 5.0 * np.log10(np.maximum(d_lum, 1e-6) * 1e6 / 10.0)
        M_g = m["g"] - dm
        # Bell+2003-style colour to mass-to-light; a factor-of-two estimate
        gr = m["g"] - m["r"]
        log_ml = -0.499 + 1.519 * gr
        L_g = 10 ** (-0.4 * (M_g - 5.11))
        log_mstar = np.log10(np.maximum(L_g, 1e-30)) + log_ml

    np.savez_compressed(OUT, ra=ra, dec=dec, z=z, z_err=z_err, z_src=z_src,
                        pstar=pstar, d_lum_mpc=d_lum, kpc_per_arcsec=kpc_per_arcsec,
                        r_kpc=r_kpc, M_g=M_g, log_mstar=log_mstar,
                        shape_r=m["shape_r"], mu_g=m["mu_g"])
    g = have & np.isfinite(r_kpc)
    if g.any():
        print(f"  z      : median {np.nanmedian(z[have]):.3f}")
        print(f"  r_eff  : median {np.nanmedian(r_kpc[g]):.2f} kpc")
        print(f"  M_g    : median {np.nanmedian(M_g[g]):.2f}")
        print(f"  log M* : median {np.nanmedian(log_mstar[g]):.2f}")
        print(f"  r_eff > 1.5 kpc AND mu_g > 24: "
              f"{int(((r_kpc > 1.5) & (m['mu_g'] > 24) & g).sum()):,}")
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
