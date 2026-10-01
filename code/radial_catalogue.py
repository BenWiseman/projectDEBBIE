"""Extract the radial representation for every catalogue object, resumably.

## Why not just run radial_stage2.py

That script builds a fixed worklist, holds every profile in a dict and writes
one npz at the end. Over 24,049 objects that is fine. Over 161,721 it is about
15 GB held in RAM for a job that runs for days, and a crash at hour 30 loses
all of it. Every long job in this project has died at least once.

So: fixed deterministic order, fixed-size shards, one npz per shard, and a
restart skips whatever is already on disk. Interrupt it freely.

## What this buys the paper

The representation result is currently measured on the DeepShadows benchmark
and never applied to the catalogue, so the two halves of the paper do not touch
and a referee can fairly ask why the method section is there. Scoring the
catalogue with the same representation joins them, and it is what makes a
two-model agreement screen possible: boosted trees and the bidirectional GRU
are two architectures over one representation, so "dual" means something
whichever way a referee cuts it.

## Rate and politeness

cutouts.fetch returns the cached file when it exists, so the ~14% of the
catalogue already on the Ironwolf drive costs nothing. For the rest the
REQUEST_DELAY_S of 1.0 s is left exactly as it is: it is a social limit against
a free service, not a technical one. RADIAL_WORKERS defaults to the same 24 the
existing stage-2 run used; do not raise it to finish sooner.

    .venv/bin/python sweeps/radial_catalogue.py            # run, resuming
    .venv/bin/python sweeps/radial_catalogue.py --probe 60 # measure the rate only
    .venv/bin/python sweeps/radial_catalogue.py --status   # how far along
"""
from __future__ import annotations

import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

sys.path.insert(0, "/home/ben/star")
sys.path.insert(0, "/home/ben/star/src")

CAT = "/home/ben/star/catalogue_debbie/v1.2_main.npy"
CACHE = "/run/media/ben/Ironwolf/Stars/project/data/cache"
QUERY_DEG = 0.25        # the service returns nothing below ~0.15; see frames.stack_at
OUTDIR = "/home/ben/legacy/radial_cat"
BENCHDIR = "/home/ben/legacy/radial_bench"
BENCH_SRC = "/home/ben/legacy/radial_profiles2.npz"
SHARD = 2000
N_WORKERS = int(os.environ.get("RADIAL_WORKERS", "24"))
ARRAYS = ("ring", "cnt", "geo", "masked_frac", "sigma")


def cached_brick_centres():
    """Centres of bricks that already have g, r and z on the Ironwolf drive."""
    import re
    have = {}
    for f in os.listdir(CACHE):
        m = re.match(r"^(\d{4}[pm]\d{3})_([grizu])_image_", f)
        if m:
            have.setdefault(m.group(1), set()).add(m.group(2))
    full = sorted(b for b, bands in have.items() if {"g", "r", "z"} <= bands)
    ra = np.array([int(b[:4]) / 10.0 for b in full])
    dec = np.array([(1 if b[4] == "p" else -1) * int(b[5:8]) / 10.0 for b in full])
    return ra, dec


def worklist(cached_only=False):
    """Catalogue objects in a fixed order, so shards are reproducible.

    cached_only keeps those sitting inside a brick already on disk. The image
    bytes are then local; only the small SIA metadata query goes out. Brick
    membership is tested geometrically because deriving a brick NAME from a
    position needs the survey's own variable-RA grid, and a hand-rolled version
    of it round-tripped only 79 of 300 real brick names.
    """
    cat = np.load(CAT, allow_pickle=True)
    ra = np.asarray(cat["ra"], float)
    dec = np.asarray(cat["dec"], float)
    keep = np.ones(len(ra), bool)
    if cached_only:
        from scipy.spatial import cKDTree

        def xyz(a, d):
            a, d = np.radians(a), np.radians(d)
            return np.column_stack([np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)])

        bra, bdec = cached_brick_centres()
        dist, _ = cKDTree(xyz(bra, bdec)).query(xyz(ra, dec), k=1)
        sep = np.degrees(2 * np.arcsin(np.clip(dist / 2, 0, 1)))
        keep = sep <= 0.125          # inside half a 0.25 deg brick
        print(f"cached bricks {len(bra):,}; objects inside one: {keep.sum():,} "
              f"of {len(ra):,} ({100 * keep.mean():.1f}%)")
    idx = np.flatnonzero(keep)
    order = idx[np.lexsort((dec[idx], ra[idx]))]
    return [{"object_id": f"{ra[i]:.5f}{dec[i]:+.5f}", "ra": float(ra[i]),
             "dec": float(dec[i]), "label": -1} for i in order]


OUT = {"dir": OUTDIR}


def shard_path(n):
    return os.path.join(OUT["dir"], f"shard_{n:04d}.npz")


def bench_worklist():
    """The 24,049 benchmark objects, re-extracted through the SAME path.

    This runs FIRST and it is the whole reason the order matters. The extraction
    path changed when the service stopped serving 0.05 deg boxes, and a spot
    check reproduced the August profiles only to corr 0.963. So the model
    trained on the August profiles cannot be trusted against anything extracted
    now. Re-extracting the benchmark and retraining makes training and
    application consistent by construction -- and if the representation does not
    survive the new extraction, that shows up after two days instead of fifteen.
    """
    z = np.load(BENCH_SRC)
    ids, ra, dec = z["object_id"], z["ra"], z["dec"]
    order = np.lexsort((dec, ra))
    return [{"object_id": str(ids[i]), "ra": float(ra[i]), "dec": float(dec[i]),
             "label": -1} for i in order]


def done_shards(total):
    return {n for n in range((total + SHARD - 1) // SHARD)
            if os.path.exists(shard_path(n)) and os.path.getsize(shard_path(n)) > 0}


def run_shard(recs, extract_one):
    out, fail = {}, 0
    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(extract_one, r): r for r in recs}
        for fut in as_completed(futs):
            r = futs[fut]
            try:
                res = fut.result()
            except Exception:                                   # noqa: BLE001
                res = None
            if res is None:
                fail += 1
            else:
                res["ra"], res["dec"] = r["ra"], r["dec"]
                out[r["object_id"]] = res
    return out, fail


def save_shard(n, out):
    keys = sorted(out)
    payload = {"object_id": np.array(keys),
               "ra": np.array([out[k]["ra"] for k in keys]),
               "dec": np.array([out[k]["dec"] for k in keys])}
    for a in ARRAYS:
        if a in out[keys[0]]:
            payload[a] = np.stack([out[k][a] for k in keys])
    tmp = shard_path(n) + ".partial"
    np.savez_compressed(tmp, **payload)
    os.replace(tmp, shard_path(n))


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    if "--benchmark" in sys.argv:
        OUT["dir"] = BENCHDIR
        os.makedirs(BENCHDIR, exist_ok=True)
        recs = bench_worklist()
    else:
        recs = worklist(cached_only="--cached-only" in sys.argv)
    total = len(recs)
    n_shards = (total + SHARD - 1) // SHARD
    have = done_shards(total)

    if "--status" in sys.argv:
        n_obj = 0
        for n in sorted(have):
            try:
                n_obj += len(np.load(shard_path(n))["object_id"])
            except Exception:                                   # noqa: BLE001
                pass
        print(f"shards {len(have):,}/{n_shards:,}   objects extracted {n_obj:,}/{total:,} "
              f"({100 * n_obj / total:.1f}%)")
        return 0

    # Ask wide and crop back: the service returns nothing at the 0.05 deg the
    # extractor wants, and simply enlarging the frame would change the zero
    # point it takes over the whole frame.
    from starfinder import frames
    _orig = frames.stack_at
    frames.stack_at = lambda *a, **k: _orig(*a, **{**k, "query_deg": QUERY_DEG})
    from train_radial_recurrent import extract_one

    if "--probe" in sys.argv:
        k = int(sys.argv[sys.argv.index("--probe") + 1])
        sample = recs[:k]
        t0 = time.time()
        out, fail = run_shard(sample, extract_one)
        el = time.time() - t0
        rate = k / el
        print(f"probe {k} objects in {el:.1f}s -> {rate:.2f}/s  "
              f"(ok {len(out)}, failed {fail})")
        print(f"projected for {total:,}: {total / rate / 3600:.1f} h "
              f"at {N_WORKERS} workers")
        if out:
            a = out[next(iter(out))]
            print("arrays returned:", {k2: np.shape(a[k2]) for k2 in ARRAYS if k2 in a})
        return 0

    print(f"{total:,} objects, {n_shards:,} shards of {SHARD}, "
          f"{len(have):,} already done, {N_WORKERS} workers", flush=True)
    t0, started = time.time(), 0
    for n in range(n_shards):
        if n in have:
            continue
        chunk = recs[n * SHARD:(n + 1) * SHARD]
        ts = time.time()
        out, fail = run_shard(chunk, extract_one)
        if out:
            save_shard(n, out)
        started += 1
        el = time.time() - t0
        left = n_shards - len(have) - started
        print(f"shard {n:04d}  ok={len(out):,} fail={fail:,}  "
              f"{len(chunk) / (time.time() - ts):.2f}/s  "
              f"elapsed {el / 3600:.1f}h  {left} shards left  "
              f"eta {left * (el / max(started, 1)) / 3600:.1f}h", flush=True)
    print("complete", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
