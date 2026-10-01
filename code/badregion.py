"""Flag candidates that come from a bad imaging region rather than the sky.

Three of twenty-four recovered candidates were the same solid-blue bad-data
square at RA 144.80, Dec -2.79, within about one arcmin of each other. A single
bad brick manufactures a whole cluster of false detections, and each one looks
like an independent candidate in a list.

Real low-surface-brightness galaxies are not clustered on arcminute scales at
this surface density -- the pool averages far less than one per square arcmin.
So an over-density of candidates in a tiny area is evidence about the *imaging*,
not about the sky.

Like the shred flag, this marks for review rather than deleting: a genuine
compact group would look the same, and those are interesting in their own right.
"""
import numpy as np

RADIUS_ARCMIN = 1.5
MIN_NEIGHBOURS = 3          # this many others within RADIUS => suspect region

# Regions confirmed bad by eye. The density test alone is not enough: once a
# clump has been thinned by later filters, the few survivors no longer have
# MIN_NEIGHBOURS around them and sail through. Two objects from the RA 144.8
# Dec -2.79 bad-data square reached the final no-designation shortlist that
# way. A region confirmed bad by looking at it stays bad.
#   (ra_deg, dec_deg, radius_arcmin, note)
BLACKLIST = [
    (144.79, -2.79, 6.0, "solid blue bad-data square; 55 fake candidates"),
]


def _unit(ra, dec):
    r, d = np.radians(ra), np.radians(dec)
    return np.column_stack([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)])


def in_blacklist(ra, dec):
    """True where a candidate falls in a region confirmed bad by eye."""
    ra = np.asarray(ra, dtype=np.float64)
    dec = np.asarray(dec, dtype=np.float64)
    bad = np.zeros(ra.shape, dtype=bool)
    for bra, bdec, brad, _ in BLACKLIST:
        sep = np.hypot((ra - bra) * np.cos(np.radians(bdec)), dec - bdec) * 60.0
        bad |= sep <= brad
    return bad


def flag(ra, dec, radius_arcmin=RADIUS_ARCMIN, min_neighbours=MIN_NEIGHBOURS):
    """True where a candidate is in an over-dense clump OR a blacklisted region."""
    from scipy.spatial import cKDTree
    ra = np.asarray(ra, dtype=np.float64)
    dec = np.asarray(dec, dtype=np.float64)
    tree = cKDTree(_unit(ra, dec))
    chord = 2.0 * np.sin(np.radians(radius_arcmin / 60.0) / 2.0)
    counts = np.array(tree.query_ball_point(_unit(ra, dec), chord,
                                            return_length=True))
    return ((counts - 1) >= min_neighbours) | in_blacklist(ra, dec)


if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else "/home/ben/legacy/clean45.npy"
    c = np.load(path)
    bad = flag(c["ra"], c["dec"])
    print(f"{path}: {len(c):,} candidates")
    print(f"  in over-dense clumps (>={MIN_NEIGHBOURS} others within "
          f"{RADIUS_ARCMIN}'): {bad.sum():,} ({bad.mean()*100:.1f}%)")
    print(f"  clean: {(~bad).sum():,}")
    np.save(path.replace(".npy", "_noclump.npy"), c[~bad])
    np.save(path.replace(".npy", "_clumped.npy"), c[bad])
    if bad.any():
        cl = c[bad]
        print("\n  worst clumps (position, count within 1.5'):")
        from scipy.spatial import cKDTree
        t = cKDTree(_unit(cl["ra"], cl["dec"]))
        ch = 2.0 * np.sin(np.radians(RADIUS_ARCMIN / 60.0) / 2.0)
        n = np.array(t.query_ball_point(_unit(cl["ra"], cl["dec"]), ch,
                                        return_length=True))
        for i in np.argsort(-n)[:6]:
            print(f"    {cl['ra'][i]:9.4f} {cl['dec'][i]:+8.4f}   {n[i]} candidates")
