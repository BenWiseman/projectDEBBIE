"""Select extended low-surface-brightness galaxy candidates straight from the
DESI Legacy Survey DR10 sweep catalogues.

Why this exists
---------------
Our own detector recovers ~24% of catalogued galaxies overall and ~1% above 40
arcsec.  The sweeps are the Legacy Survey team's own measurements of the *same
imaging we have been searching*: position, morphological type, fluxes, and
`shape_r` (the half-light radius of the fitted profile, in arcsec).  Selecting
candidates from the sweep catalogue skips our detector entirely.

What this does NOT do
---------------------
It does not find anything the Legacy Survey pipeline itself failed to detect.
Tractor has its own LSB incompleteness.  This is a different, larger net --
not a complete one.  Objects selected here are real detections; whether they
are *new* is a separate question answered by the catalogue screen.

Selection
---------
UDG-like, following van Dokkum et al. 2015 in spirit:
  - type in {EXP, DEV, SER, REX}  (resolved, not PSF -- excludes stars)
  - shape_r >= 3 arcsec           (extended)
  - mu_eff(g) >= 24 mag/arcsec^2  (low surface brightness)
  - detected in g and r           (so a colour exists)

mu_eff = m + 2.5*log10(2 * pi * r_e^2), the mean SB within the half-light
radius, which is what shape_r and the total flux give us directly.
"""
import glob
import os
import sys

import numpy as np
from astropy.io import fits

RESOLVED = {"EXP", "DEV", "SER", "REX"}
MIN_R = 3.0          # arcsec
MIN_MU = 24.0        # mag/arcsec^2, mean SB inside r_e
# Above mu ~26.5 Tractor is fitting a large radius to noise: 62% of the first
# pool sat there and the contact sheets for it were almost entirely empty sky
# with a crosshair. Real objects in these data live at mu 24-26.5.
MAX_MU = 26.5

COLS = ["RA", "DEC", "TYPE", "SHAPE_R", "SHAPE_E1", "SHAPE_E2",
        "FLUX_G", "FLUX_R", "FLUX_Z", "FLUX_IVAR_G", "FLUX_IVAR_R",
        "NOBS_G", "NOBS_R", "MASKBITS", "BRICKNAME", "OBJID"]


def nanomaggie_to_mag(f):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(f > 0, 22.5 - 2.5 * np.log10(np.maximum(f, 1e-12)), np.nan)


def mine(path):
    """Return a structured array of LSB candidates from one sweep file."""
    with fits.open(path, memmap=True) as h:
        d = h[1].data
        typ = np.char.strip(d["TYPE"].astype(str))
        r = d["SHAPE_R"].astype(np.float64)

        keep = np.isin(typ, list(RESOLVED)) & (r >= MIN_R)
        # require real coverage in g and r, and no bright-star / bleed mask
        keep &= (d["NOBS_G"] > 0) & (d["NOBS_R"] > 0)
        keep &= d["MASKBITS"] == 0
        keep &= (d["FLUX_G"] > 0) & (d["FLUX_R"] > 0)
        # Single-epoch regions are where the bad-data squares and satellite
        # trails survive into the coadd: with only one exposure there is no
        # second frame to reject the outlier against. Requiring two epochs in
        # both bands removes most of the artefacts seen on the contact sheets.
        keep &= (d["NOBS_G"] >= 2) & (d["NOBS_R"] >= 2)
        # A real galaxy is detected in z as well; trails and bad columns are
        # usually confined to one or two bands.
        if "FLUX_Z" in d.columns.names:
            keep &= d["FLUX_Z"] > 0
        if not keep.any():
            return None

        idx = np.flatnonzero(keep)
        g = nanomaggie_to_mag(d["FLUX_G"][idx])
        rm = nanomaggie_to_mag(d["FLUX_R"][idx])
        re = r[idx]
        # mean surface brightness within the half-light radius
        mu_g = g + 2.5 * np.log10(2.0 * np.pi * re ** 2)

        sel = (mu_g >= MIN_MU) & (mu_g <= MAX_MU) & np.isfinite(mu_g)
        if not sel.any():
            return None
        idx = idx[sel]

        out = np.zeros(idx.size, dtype=[
            ("ra", "f8"), ("dec", "f8"), ("shape_r", "f4"),
            ("g", "f4"), ("r", "f4"), ("mu_g", "f4"),
            ("type", "U4"), ("brick", "U16"), ("objid", "i8")])
        out["ra"] = d["RA"][idx]
        out["dec"] = d["DEC"][idx]
        out["shape_r"] = re[sel]
        out["g"] = g[sel]
        out["r"] = rm[sel]
        out["mu_g"] = mu_g[sel]
        out["type"] = typ[idx]
        out["brick"] = np.char.strip(d["BRICKNAME"][idx].astype(str))
        out["objid"] = d["OBJID"][idx]
        return out


def complete(path):
    """A sweep is only usable once fully downloaded."""
    try:
        with fits.open(path) as h:
            h[1].data[:1]
        return True
    except Exception:
        return False


if __name__ == "__main__":
    files = sorted(glob.glob("/home/ben/legacy/sweeps/*.fits"))
    if len(sys.argv) > 1:
        files = files[: int(sys.argv[1])]
    chunks, n_files, n_rows = [], 0, 0
    for f in files:
        if not complete(f):
            continue
        try:
            c = mine(f)
        except Exception as e:
            print(f"  {os.path.basename(f)}: FAILED {type(e).__name__}: {e}")
            continue
        n_files += 1
        with fits.open(f, memmap=True) as h:
            n_rows += len(h[1].data)
        if c is not None:
            chunks.append(c)
            print(f"  {os.path.basename(f)}: {len(c)} LSB candidates")
    if not chunks:
        print("no complete sweep files yet")
        raise SystemExit(0)
    allc = np.concatenate(chunks)
    np.save("/home/ben/legacy/lsb_candidates.npy", allc)
    area = n_files * 25.0  # 5x5 deg boxes, ignoring survey edges
    print(f"\n{n_files} sweeps, {n_rows:,} sources -> {len(allc):,} LSB candidates")
    print(f"~{len(allc)/max(area,1):.1f} per deg^2 over ~{area:.0f} deg^2 (box area, not masked area)")
    print(f"shape_r: median {np.median(allc['shape_r']):.1f}\" "
          f"90th {np.percentile(allc['shape_r'],90):.1f}\" max {allc['shape_r'].max():.1f}\"")
    print(f"mu_g   : median {np.median(allc['mu_g']):.1f}")
    print("saved /home/ben/legacy/lsb_candidates.npy")
