"""Distil every sweep into one compact table, so the raw sweeps become optional.

The sweeps are 2.44 TB and both disks are full. But 99% of that is sources we
will never reconsider: point sources, and resolved objects far below any size
cut we would plausibly use. Keeping every resolved source with r_eff >= 2
arcsec -- comfortably under the current 3 arcsec threshold -- is 55 million rows
in about 2.6 GB. That is 924x smaller and preserves the thing the sweeps were
being kept for.

Why the sweeps were kept at all
-------------------------------
The mining criteria changed five times in one session (MAX_MU 28.5 -> 26.5,
NOBS >= 2, FLUX_Z > 0, shred radius, bad-region filter). Deleting the source
data would have destroyed the ability to re-apply our own corrections. This
extract keeps that ability: any future change to the surface-brightness cut,
colour cut, morphology selection or quality flags can be re-run against it in
seconds, with no download.

What it does NOT preserve
-------------------------
- Sources with r_eff < 2 arcsec, or PSF-type sources. If a future selection
  wants those, the sweeps must be re-downloaded.
- Columns not listed in KEEP below.
Both limits are deliberate and stated so nobody later assumes this is lossless.

Nothing is deleted by this script. It only writes the extract; removing or
moving sweeps is a separate, explicit decision.
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np
from astropy.io import fits

LEGACY = "/home/ben/legacy"
OUT = os.path.join(LEGACY, "resolved_sources.npz")
SWEEP_DIRS = [f"{LEGACY}/sweeps", "/run/media/ben/Fanxiang/sweeps"]
MIN_R = 2.0          # below the 3.0 arcsec production cut, with headroom

# Columns are kept as raw fluxes, not magnitudes, so nothing is lost to a
# conversion we might later want to redo differently.
#
# The first version of this kept 15 of 130 columns and would have silently
# capped every future search: no flux uncertainties (hence no error bars on any
# magnitude), no EBV (hence no galactic extinction correction), no FRACFLUX
# (hence no blending cut, and blending is a failure mode we already hit), and no
# BRICKID/OBJID (hence only positional joins to photo-z instead of exact ones).
# Widening it roughly triples the size, which against 2.44 TB of sweeps is
# nothing.
DTYPE = [("ra", "f8"), ("dec", "f8"),
         ("brickid", "i4"), ("objid", "i4"),          # exact join key
         ("shape_r", "f4"), ("shape_r_ivar", "f4"),
         ("e1", "f4"), ("e2", "f4"), ("sersic", "f4"),
         ("flux_g", "f4"), ("flux_r", "f4"), ("flux_z", "f4"),
         ("flux_w1", "f4"), ("flux_w2", "f4"),
         ("ivar_g", "f4"), ("ivar_r", "f4"), ("ivar_z", "f4"),
         ("ivar_w1", "f4"),                            # error bars
         ("fracflux_g", "f4"), ("fracflux_r", "f4"),   # blending
         ("fracmasked_g", "f4"), ("fracmasked_r", "f4"),
         ("fracin_g", "f4"),
         ("ebv", "f4"),                                # extinction
         ("dchisq_psf", "f4"), ("dchisq_ser", "f4"),   # morphology preference
         ("nobs_g", "i2"), ("nobs_r", "i2"), ("nobs_z", "i2"),
         ("maskbits", "i4"), ("type", "U4")]


def _col(d, name, idx, default=np.nan):
    """Read a column if the release has it, else fill. DR10 sweeps vary."""
    if name in d.columns.names:
        return d[name][idx]
    return np.full(len(idx), default)


def mag(f):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(f > 0, 22.5 - 2.5 * np.log10(np.maximum(f, 1e-12)), np.nan)


def main():
    files = sorted(f for d in SWEEP_DIRS for f in glob.glob(f"{d}/*.fits"))
    print(f"{len(files)} sweeps to distil (r_eff >= {MIN_R}\")", flush=True)
    chunks, n_src, n_bad = [], 0, 0

    for k, f in enumerate(files, 1):
        try:
            with fits.open(f, memmap=True) as h:
                d = h[1].data
                n_src += len(d)
                typ = np.char.strip(d["TYPE"].astype(str))
                keep = (typ != "PSF") & (d["SHAPE_R"] >= MIN_R)
                idx = np.flatnonzero(keep)
                if not idx.size:
                    continue
                out = np.zeros(idx.size, dtype=DTYPE)
                out["ra"] = d["RA"][idx]
                out["dec"] = d["DEC"][idx]
                out["brickid"] = _col(d, "BRICKID", idx, 0)
                out["objid"] = _col(d, "OBJID", idx, 0)
                out["shape_r"] = d["SHAPE_R"][idx]
                out["shape_r_ivar"] = _col(d, "SHAPE_R_IVAR", idx)
                out["e1"] = d["SHAPE_E1"][idx]
                out["e2"] = d["SHAPE_E2"][idx]
                out["sersic"] = _col(d, "SERSIC", idx)
                for b in ("G", "R", "Z", "W1", "W2"):
                    out[f"flux_{b.lower()}"] = _col(d, f"FLUX_{b}", idx)
                for b in ("G", "R", "Z", "W1"):
                    out[f"ivar_{b.lower()}"] = _col(d, f"FLUX_IVAR_{b}", idx)
                for b in ("G", "R"):
                    out[f"fracflux_{b.lower()}"] = _col(d, f"FRACFLUX_{b}", idx)
                    out[f"fracmasked_{b.lower()}"] = _col(d, f"FRACMASKED_{b}", idx)
                out["fracin_g"] = _col(d, "FRACIN_G", idx)
                out["ebv"] = _col(d, "EBV", idx)
                dq = _col(d, "DCHISQ", idx, np.nan)
                if getattr(dq, "ndim", 1) == 2 and dq.shape[1] >= 5:
                    out["dchisq_psf"] = dq[:, 0]
                    out["dchisq_ser"] = dq[:, 4]
                out["nobs_g"] = d["NOBS_G"][idx]
                out["nobs_r"] = d["NOBS_R"][idx]
                out["nobs_z"] = _col(d, "NOBS_Z", idx, 0)
                out["maskbits"] = d["MASKBITS"][idx]
                out["type"] = typ[idx]
                chunks.append(out)
        except Exception as e:                                    # noqa: BLE001
            n_bad += 1
            print(f"  {os.path.basename(f)}: {type(e).__name__}", flush=True)
            continue
        if k % 100 == 0:
            have = sum(len(c) for c in chunks)
            print(f"  {k}/{len(files)}  {have:,} rows kept", flush=True)

    if not chunks:
        print("nothing extracted")
        return
    allr = np.concatenate(chunks)
    np.savez_compressed(OUT, rows=allr, min_r=MIN_R,
                        n_sweeps=len(files) - n_bad, n_sources_scanned=n_src)
    size = os.path.getsize(OUT) / 1e9
    print(f"\n{len(allr):,} resolved sources from {n_src:,} scanned")
    print(f"  {n_bad} sweeps unreadable (truncated)")
    print(f"  -> {OUT}  {size:.2f} GB")
    print(f"  sweeps on disk are 2.44 TB, i.e. {2440/max(size,0.01):.0f}x larger")


if __name__ == "__main__":
    main()
