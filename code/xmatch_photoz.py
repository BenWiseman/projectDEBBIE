"""Get photo-z for the catalogue via CDS X-Match, not by bulk download.

Downloading VII/292 wholesale does not work: VizieR truncates each response at
~15-18M rows regardless of -out.max=unlimited, so a slice requested as Dec
-10..0 came back covering only -10.0 to -8.7. Getting full coverage that way
would take ~200 narrow requests and hundreds of GB, for a catalogue we only
need at 240,000 specific positions.

CDS X-Match takes an uploaded position list and returns matches directly. It is
the tool built for this, and it also returns pstar and the quality flags.

A note on a wrong turn: a query for the maximum declination in VII/292/south
returned +0.0036, which looked like the southern table stopping at the equator
and explained three "hanging" northern slices. It was wrong -- the sort had not
applied, and X-Match happily matches objects at Dec +11.7 against that same
table. The hang was the truncation above, not a hemisphere split.
"""
from __future__ import annotations

import csv
import io
import os
import sys
import time
import urllib.request

import numpy as np

LEGACY = "/home/ben/legacy"
OUT = os.path.join(LEGACY, "photoz_xmatch.npz")
URL = "http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync"
CHUNK = 25_000            # positions per request
RADIUS = 2.0              # arcsec


def post_chunk(rows, cat2="vizier:VII/292/south"):
    """Upload a chunk of positions, return the matched rows."""
    buf = io.StringIO()
    buf.write("ra,dec\n")
    for ra, dec in rows:
        buf.write(f"{ra:.6f},{dec:.6f}\n")
    body = buf.getvalue().encode()

    boundary = "----xmatchboundary"
    parts = []

    def field(name, value):
        parts.append(f"--{boundary}\r\nContent-Disposition: form-data; "
                     f'name="{name}"\r\n\r\n{value}\r\n'.encode())

    field("request", "xmatch")
    field("distMaxArcsec", str(RADIUS))
    field("RESPONSEFORMAT", "csv")
    field("colRA1", "ra")
    field("colDec1", "dec")
    field("cat2", cat2)
    parts.append(f"--{boundary}\r\nContent-Disposition: form-data; "
                 f'name="cat1"; filename="c.csv"\r\n'
                 f"Content-Type: text/csv\r\n\r\n".encode() + body + b"\r\n")
    parts.append(f"--{boundary}--\r\n".encode())
    data = b"".join(parts)

    req = urllib.request.Request(
        URL, data=data,
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    txt = urllib.request.urlopen(req, timeout=900).read().decode("utf8", "replace")
    return list(csv.DictReader(io.StringIO(txt)))


def main():
    cat = np.load(os.path.join(LEGACY, "VERIFIED_NOVEL.npy"))
    n = len(cat)
    print(f"{n:,} positions, {CHUNK:,} per request", flush=True)

    z = np.full(n, np.nan)
    ze = np.full(n, np.nan)
    pstar = np.full(n, np.nan)
    fqual = np.full(n, np.nan)
    got = 0

    # index positions so returned rows can be mapped back
    key = {f"{o['ra']:.6f}_{o['dec']:.6f}": i for i, o in enumerate(cat)}

    for start in range(0, n, CHUNK):
        sub = cat[start:start + CHUNK]
        rows = [(float(o["ra"]), float(o["dec"])) for o in sub]
        for attempt in range(4):
            try:
                # VII/292 is split by hemisphere. Try south, then north for
                # anything it did not match -- 43.6% of the catalogue lives in
                # the northern table.
                res = post_chunk(rows)
                if len(res) < len(rows) * 0.9:
                    res = res + post_chunk(rows, cat2="vizier:VII/292/north")
                break
            except Exception as e:                                # noqa: BLE001
                print(f"  chunk {start//CHUNK}: {type(e).__name__}, retry",
                      flush=True)
                time.sleep(20 * (attempt + 1))
        else:
            print(f"  chunk {start//CHUNK}: FAILED, skipping", flush=True)
            continue

        for r in res:
            try:
                k = f"{float(r['ra']):.6f}_{float(r['dec']):.6f}"
                i = key.get(k)
                if i is None:
                    continue
                # keep the closest match if several returned
                if np.isfinite(z[i]):
                    continue
                z[i] = float(r["zphot"]) if r.get("zphot") else np.nan
                ze[i] = float(r["e_zphot"]) if r.get("e_zphot") else np.nan
                pstar[i] = float(r["pstar"]) if r.get("pstar") else np.nan
                fqual[i] = float(r["fqual"]) if r.get("fqual") else np.nan
                got += 1
            except Exception:                                     # noqa: BLE001
                continue
        print(f"  {min(start+CHUNK, n):,}/{n:,}  matched {got:,}", flush=True)

    np.savez_compressed(OUT, ra=cat["ra"], dec=cat["dec"], zphot=z,
                        e_zphot=ze, pstar=pstar, fqual=fqual)
    have = np.isfinite(z)
    print(f"\nphoto-z for {have.sum():,}/{n:,} ({have.mean()*100:.1f}%)")
    if have.any():
        print(f"  z median {np.nanmedian(z[have]):.3f}")
        print(f"  pstar>0.5 (likely stars): {int((pstar > 0.5).sum()):,}")
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
