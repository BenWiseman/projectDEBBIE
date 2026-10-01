"""Export the catalogue as a real, publishable artifact.

A .npy file is not a catalogue. A catalogue is: a designation per object, a
position anyone can point a telescope at, the measurements, the provenance of
every check that was run, and a machine-readable table in a format journals and
CDS accept. This writes that.

Designations follow ACRONYM (below): `WINNIE NNN`, sequential, no
zero padding (NNN is a maximum width, not a fixed one). Numbers are assigned in
order of first appearance and are STABLE -- once an object has a number it keeps
it, because a designation that moves between versions is worse than useless.

Outputs
-------
    catalogue/<ACRONYM>.csv  the table (CDS-friendly, one row per object)
    catalogue/<ACRONYM>.ecsv same with full column metadata and units
    catalogue/README.md      what the columns mean and how objects were selected
    catalogue/assignments.json  designation -> position, append-only and stable
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, "/home/ben/star/src")
sys.path.insert(0, "/home/ben/legacy")

LEGACY = "/home/ben/legacy"
# The designation prefix. NOTE: an acronym must be registered with the CDS
# Dictionary of Nomenclature before it can be used in a publication. DEBBIE is
# registered; DEBBIE is the one to file.
# DEBBIE = Deep Extragalactic Broad-Band Inference Experiment. Ben's decision
# 2026-09-04: DEBBIE is the acronym that will be registered with the CDS
# Dictionary of Nomenclature, so the released designations carry it. Numbers
# are keyed by position in assignments.json and never change across renames
# (both earlier acronyms preserved them).
ACRONYM = "DEBBIE"
OUT = "/home/ben/star/catalogue_winnie"
ASSIGN = os.path.join(OUT, "assignments.json")


def sexagesimal(ra, dec):
    """CDS position format: 'hh mm ss.ss +dd mm ss.s'."""
    h = ra / 15.0
    hh = int(h); mm = int((h - hh) * 60); ss = (h - hh - mm / 60.0) * 3600
    sign = "+" if dec >= 0 else "-"
    a = abs(dec)
    dd = int(a); am = int((a - dd) * 60); asec = (a - dd - am / 60.0) * 3600
    return f"{hh:02d} {mm:02d} {ss:05.2f} {sign}{dd:02d} {am:02d} {asec:04.1f}"


def iau_name(ra, dec):
    """The positional form, for cross-referencing: <ACRONYM> Jhhmmss.s+ddmmss."""
    h = ra / 15.0
    hh = int(h); mm = int((h - hh) * 60); ss = (h - hh - mm / 60.0) * 3600
    sign = "+" if dec >= 0 else "-"
    a = abs(dec)
    dd = int(a); am = int((a - dd) * 60); asec = (a - dd - am / 60.0) * 3600
    return f"{ACRONYM} J{hh:02d}{mm:02d}{ss:04.1f}{sign}{dd:02d}{am:02d}{asec:02.0f}"


def load_assignments():
    if os.path.exists(ASSIGN):
        try:
            return json.load(open(ASSIGN))
        except json.JSONDecodeError:
            pass
    return {}


def key(ra, dec):
    return f"{ra:.5f}_{dec:+.5f}"


def rebuild_verified_novel():
    """Regenerate VERIFIED_NOVEL.npy from the live verification state.

    This used to be a hand-run step, so the catalogue froze at 6,898 while
    verification had reached 9,752 novel. Unattended for a week that means the
    catalogue simply stops growing while the pipeline keeps working -- the one
    artifact that has to exist to publish, stuck.

    Re-screening here also catches objects verified earlier that a since-
    tightened filter now rejects, so the count reflects current criteria rather
    than whatever was true when each object happened to be checked.
    """
    import badregion
    from starfinder import localcat

    state = os.path.join(LEGACY, "verify_state.json")
    if not os.path.exists(state):
        return None
    st = json.load(open(state))
    rows = [(float(k.split("_")[0]), float(k.split("_")[1]), v["r"], v["mu"])
            for k, v in st.items() if v["v"] == "novel"]
    if not rows:
        return None
    a = np.array(rows)
    cat = localcat.load()
    ok = (~cat.known(a[:, 0], a[:, 1], 20.0)
          & ~cat.shredded(a[:, 0], a[:, 1], 45.0)
          & ~badregion.flag(a[:, 0], a[:, 1]))
    keep = a[ok]
    out = np.zeros(len(keep), dtype=[("ra", "f8"), ("dec", "f8"),
                                     ("shape_r", "f4"), ("mu_g", "f4"),
                                     ("type", "U8")])
    out["ra"], out["dec"] = keep[:, 0], keep[:, 1]
    out["shape_r"], out["mu_g"] = keep[:, 2], keep[:, 3]
    out = out[np.argsort(-out["shape_r"])]

    # Collapse entries that are the same galaxy.
    #
    # Tractor shreds a galaxy into several sources, and the shred filter only
    # compares against NAMED galaxies -- never against other catalogue entries.
    # So one galaxy could receive several WINNIE designations: 172.7337-13.6117
    # and 172.7330-13.6134 are one object 7 arcsec apart, both in the figure
    # shortlist. 2.2% of entries had a neighbour within 10 arcsec.
    #
    # The radius scales with size: two entries 7 arcsec apart on a 13 arcsec
    # object are certainly one galaxy, while two compact 3 arcsec objects that
    # far apart may genuinely be a pair. Sorted largest-first above, so the
    # survivor is the largest member -- the one most likely to be the galaxy
    # rather than a fragment of it.
    from scipy.spatial import cKDTree
    from starfinder.localcat import _unit
    tree = cKDTree(_unit(out["ra"], out["dec"]))
    keep = np.ones(len(out), dtype=bool)
    for i in range(len(out)):
        if not keep[i]:
            continue
        rad = max(6.0, float(out["shape_r"][i]))
        chord = 2.0 * np.sin(np.radians(rad / 3600.0) / 2.0)
        for j in tree.query_ball_point(_unit(out["ra"][i:i+1],
                                             out["dec"][i:i+1])[0], chord):
            if j != i:
                keep[j] = False
    n_dup = int((~keep).sum())
    out = out[keep]
    np.save(os.path.join(LEGACY, "VERIFIED_NOVEL.npy"), out)
    print(f"  re-screened {len(a):,} novel -> {len(out):,} pass current filters "
          f"({n_dup:,} merged as shreds of the same galaxy)")
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    cat = rebuild_verified_novel()
    if cat is None:
        cat = np.load(os.path.join(LEGACY, "VERIFIED_NOVEL.npy"))
    # brightest/largest first is a stable, defensible ordering for assignment
    cat = cat[np.argsort(-cat["shape_r"])]

    assigned = load_assignments()
    next_n = (max((v["n"] for v in assigned.values()), default=0) + 1)
    verify = {}
    try:
        verify = json.load(open(os.path.join(LEGACY, "verify_state.json")))
    except Exception:                                             # noqa: BLE001
        pass

    rows = []
    for o in cat:
        ra, dec = float(o["ra"]), float(o["dec"])
        k = key(ra, dec)
        if k not in assigned:
            assigned[k] = {"n": next_n, "ra": ra, "dec": dec}
            next_n += 1
        n = assigned[k]["n"]
        v = verify.get(f"{ra:.5f}_{dec:.5f}", {})
        desig = (v.get("d") or [])
        rows.append({
            "name": f"{ACRONYM} {n}",
            "iau": iau_name(ra, dec),
            "ra_deg": round(ra, 6),
            "dec_deg": round(dec, 6),
            "position": sexagesimal(ra, dec),
            "r_eff_arcsec": round(float(o["shape_r"]), 2),
            "mu_eff_g": round(float(o["mu_g"]), 2),
            "prior_designations": ";".join(str(x) for x in desig[:3]),
        })

    json.dump(assigned, open(ASSIGN, "w"), indent=0)

    # CSV
    import csv
    with open(os.path.join(OUT, f"{ACRONYM}.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ECSV with units, which is what most journals and CDS actually want
    try:
        from astropy.table import Table
        import astropy.units as u
        t = Table(rows={k: [r[k] for r in rows] for k in rows[0]}) \
            if False else Table([[r[k] for r in rows] for k in rows[0]],
                                names=list(rows[0].keys()))
        t["ra_deg"].unit = u.deg
        t["dec_deg"].unit = u.deg
        t["r_eff_arcsec"].unit = u.arcsec
        t["mu_eff_g"].unit = u.mag / u.arcsec ** 2
        t.meta["acronym"] = ACRONYM
        t.meta["survey"] = "DESI Legacy Imaging Survey DR10 (south)"
        t.write(os.path.join(OUT, f"{ACRONYM}.ecsv"), format="ascii.ecsv",
                overwrite=True)
    except Exception as e:                                        # noqa: BLE001
        print(f"  ECSV write failed: {type(e).__name__}: {e}")

    print(f"{len(rows):,} objects -> {OUT}/{ACRONYM}.csv")
    print(f"  designations {ACRONYM} 1 .. {ACRONYM} {next_n - 1}")
    print(f"  r_eff {cat['shape_r'].min():.1f}-{cat['shape_r'].max():.1f} arcsec")
    return rows


if __name__ == "__main__":
    main()
