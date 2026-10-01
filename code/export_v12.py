"""Write the v1.2 release tables the paper says it releases.

This existed only as .npy plus a numbers file: the Data Availability section
promised DEBBIE.csv at 161,721 rows and DEBBIE_lowb_extension.csv at 78,713,
and neither file existed anywhere. For a catalogue paper the catalogue IS the
result, so this is the one file that cannot be missing.

Designations come from assignments.json and are stable -- an object keeps the
number it was first given, across two earlier acronyms and across
v1.0 -> v1.2. Removed objects leave gaps; that is correct, because a
designation that gets reused points at two different objects in two papers.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, "/home/ben/star/src")
sys.path.insert(0, "/home/ben/legacy")
from starfinder.localcat import _unit                             # noqa: E402
from export_catalogue import sexagesimal, iau_name, key, ACRONYM  # noqa: E402

OUT = "/home/ben/star/catalogue_debbie/v1.2"
ASSIGN = "/home/ben/star/catalogue_winnie/assignments.json"


def build(npy, stem, assigned, verify, phys, tree):
    cat = np.load(npy)
    d, j = tree.query(_unit(cat["ra"], cat["dec"]), k=1,
                      distance_upper_bound=2 * np.sin(np.radians(1 / 3600) / 2))
    hit = np.isfinite(d)
    # Extract the npz columns ONCE. Indexing phys["z"] inside the row loop
    # re-decompresses the whole column every iteration: measured at 161s per
    # 20,000 rows, i.e. ~22 minutes for this table, versus 5 seconds hoisted.
    pz, psrc, pps = phys["z"], phys["z_src"], phys["pstar"]
    rows = []
    import time as _t; _t0=_t.time()
    for i, o in enumerate(cat):
        if i and i % 20000 == 0:
            print(f'    row {i:,} ({_t.time()-_t0:.0f}s)', flush=True)
        ra, dec = float(o["ra"]), float(o["dec"])
        a = assigned.get(key(ra, dec))
        if a is None:
            sys.exit(f"UNASSIGNED DESIGNATION for {ra} {dec}")
        v = verify.get(f"{ra:.5f}_{dec:.5f}", {})
        # NaN, not "", for absent numbers. Mixing floats and empty strings in a
        # column makes astropy fall back to an object dtype, and writing an
        # object column to ECSV is pathologically slow -- this is what turned a
        # two-second export into a ten-minute hang.
        z = float("nan"); ps = float("nan"); zs = ""
        if hit[i]:
            zz = pz[j[i]]
            if np.isfinite(zz) and zz > 0:
                z = round(float(zz), 4)
                zs = str(psrc[j[i]])
            pp = pps[j[i]]
            if np.isfinite(pp):
                ps = round(float(pp), 3)
        rows.append({
            "name": f"{ACRONYM} {a['n']}",
            "iau": iau_name(ra, dec),
            "ra_deg": round(ra, 6), "dec_deg": round(dec, 6),
            "position": sexagesimal(ra, dec),
            "r_eff_arcsec": round(float(o["shape_r"]), 2),
            "mu_eff_g": round(float(o["mu_g"]), 2),
            "z_phot": z, "z_source": zs, "pstar": ps,
            "prior_designations": ";".join(str(x) for x in (v.get("d") or [])[:3]),
        })
    print(f"    loop done ({_t.time()-_t0:.0f}s), writing csv", flush=True)
    with open(f"{OUT}/{stem}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    print("    csv written, building table", flush=True)
    from astropy.table import Table
    import astropy.units as u
    cols = {}
    for k in rows[0]:
        v = [r[k] for r in rows]
        cols[k] = np.array(v, dtype=float) if k in ("z_phot", "pstar") else v
    t = Table(list(cols.values()), names=list(cols.keys()))
    t["ra_deg"].unit = u.deg; t["dec_deg"].unit = u.deg
    t["r_eff_arcsec"].unit = u.arcsec
    t["mu_eff_g"].unit = u.mag / u.arcsec ** 2
    t.meta["acronym"] = ACRONYM
    t.meta["version"] = "1.2"
    t.meta["survey"] = "DESI Legacy Imaging Survey DR10 (south)"
    print("    writing ecsv", flush=True)
    t.write(f"{OUT}/{stem}.ecsv", format="ascii.ecsv", overwrite=True)
    return len(rows)


def main():
    print("start", flush=True)
    assigned = json.load(open(ASSIGN))
    print("assignments loaded", flush=True)
    verify = json.load(open("/home/ben/legacy/verify_state.json"))
    print("verify loaded", flush=True)
    phys = np.load("/home/ben/legacy/physical.npz", allow_pickle=True)
    print("phys loaded", flush=True)
    tree = cKDTree(_unit(phys["ra"], phys["dec"]))
    print("tree built", flush=True)
    print("building main", flush=True)
    n1 = build(f"{OUT}/../v1.2_main.npy", "DEBBIE", assigned, verify, phys, tree)
    print(f"main done: {n1}", flush=True)
    n2 = build(f"{OUT}/../v1.2_ext.npy", "DEBBIE_lowb_extension",
               assigned, verify, phys, tree)
    print(f"DEBBIE.csv {n1:,} rows;  DEBBIE_lowb_extension.csv {n2:,} rows")

    lines = []
    for f in sorted(os.listdir(OUT)):
        if f == "SHA256":
            continue
        h = hashlib.sha256(open(f"{OUT}/{f}", "rb").read()).hexdigest()
        lines.append(f"{h}  {f}")
    open(f"{OUT}/SHA256", "w").write("\n".join(lines) + "\n")
    print(f"checksummed {len(lines)} files")
    return n1, n2


if __name__ == "__main__":
    main()
