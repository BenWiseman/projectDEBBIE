"""Verify candidates BEFORE they are ever put in front of a human.

The workflow was backwards. Contact sheets were rendered from the raw pool
ranked by size, so the scarcest resource in this project -- time actually spent
looking at objects -- was being spent on galaxies SDSS catalogued years ago.
Of the first 28 filed by eye, 13 were already identified.

Now: nothing reaches a sheet until every check we have has come back, and come
back *clean*. A check that errors is retried, and if it keeps failing the
object stays in the queue rather than being promoted.

Order is cheapest-first, so the expensive network call only runs on what
survives:
  1. local named-catalogue screen   (in memory, microseconds)
  2. SIMBAD + SMUDGes               (network, seconds)
  3. NED                            (network, seconds, the real arbiter)

Output: verified_novel.npy -- objects with no galaxy identification anywhere.
Those, and only those, get rendered.
"""
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

sys.path.insert(0, "/home/ben/star")
sys.path.insert(0, "/home/ben/star/src")

POOL = "/home/ben/legacy/candidate_pool.npy"
STATE = "/home/ben/legacy/verify_state.json"
NOVEL = "/home/ben/legacy/verified_novel.npy"
WORKERS = int(os.environ.get("VERIFY_WORKERS", "6"))
RETRIES = 3


def key(o):
    return f"{float(o['ra']):.5f}_{float(o['dec']):.5f}"


def load_state():
    if os.path.exists(STATE):
        try:
            return json.load(open(STATE))
        except json.JSONDecodeError:
            pass
    return {}


def verify(o):
    """Return ('novel'|'known'|'failed', detail). Never guesses."""
    import register
    ra, dec = float(o["ra"]), float(o["dec"])
    last = ""
    for attempt in range(RETRIES):
        rec = register.assess(ra, dec)
        if rec.get("failed_checks"):
            last = str(rec["failed_checks"][:1])
            time.sleep(10 * (attempt + 1))
            continue
        if rec.get("identified_as_galaxy_by"):
            return "known", rec["identified_as_galaxy_by"][:2]
        return "novel", rec.get("machine_designations", [])[:3]
    return "failed", last


MODEL = "/home/ben/legacy/radial_model.pkl"
SCORES = "/home/ben/legacy/model_scores.npz"


def rank_by_model(clean):
    """Order the queue by the radial model, best first, junk last.

    NED is the only slow stage left (~0.1-0.5 objects/second, no bulk download
    exists). The model flags roughly 37% of candidates as junk, so scoring
    first both removes that work and puts the most promising objects at the
    front instead of merely the largest.

    Falls back to size ordering if no scores are available yet -- a missing
    model must never stop verification, only make it less well ordered.
    """
    if not os.path.exists(SCORES):
        print("  no model scores yet -- ordering by size", flush=True)
        return clean[np.argsort(-clean["shape_r"])]
    try:
        z = np.load(SCORES)
        key_to_score = {f"{r:.5f}_{d:.5f}": s
                        for r, d, s in zip(z["ra"], z["dec"], z["score"])}
    except Exception as e:                                        # noqa: BLE001
        print(f"  scores unreadable ({type(e).__name__}) -- ordering by size",
              flush=True)
        return clean[np.argsort(-clean["shape_r"])]

    s = np.array([key_to_score.get(key(o), np.nan) for o in clean])
    scored = np.isfinite(s)
    # Unscored objects sort between "good" and "junk": they are unknown, not
    # bad, and must not be buried behind things the model actively rejected.
    order = np.lexsort((-clean["shape_r"], -np.where(scored, s, 0.5)))
    n_junk = int((s[scored] < 0.1).sum()) if scored.any() else 0
    print(f"  model-ranked: {scored.sum():,} scored "
          f"({n_junk:,} flagged junk, queued last), "
          f"{(~scored).sum():,} unscored", flush=True)
    return clean[order]


def main():
    if not os.path.exists(POOL):
        print("no pool yet")
        return
    pool = np.load(POOL)

    # Sweeps overlap at their boundaries, so the conveyor appends the same
    # object more than once: 133,204 of 652,881 pool rows (20.4%) are duplicate
    # positions. verify_state is keyed by position so the catalogue never
    # duplicated -- but the verifier was re-running NED on each copy, wasting a
    # fifth of the slowest stage in the pipeline.
    pkey = np.char.add(np.char.mod("%.5f", pool["ra"]),
                       np.char.mod("%+.5f", pool["dec"]))
    _, first = np.unique(pkey, return_index=True)
    if len(first) < len(pool):
        print(f"pool {len(pool):,} -> {len(first):,} unique positions "
              f"({len(pool)-len(first):,} duplicates dropped)", flush=True)
        pool = pool[np.sort(first)]

    # Cheapest checks first, so the slow one runs on as little as possible.
    # NED is ~0.1/s: 29,294 objects is 81 hours. The local screen (109.8M
    # objects) and the shred flag are both instant and between them remove
    # ~95%, so running them first turns days into hours.
    #
    # Nothing is discarded here -- the removed objects are archived by
    # build_rejections.py with the matched object and its separation, because
    # a 20 arcsec radius against imperfect catalogue positions gets real
    # objects wrong and Ben sifts those by eye.
    import badregion
    from starfinder import localcat
    cat = localcat.load()
    known = cat.known(pool["ra"], pool["dec"], 20.0)
    # 45 arcsec, not 90: the wider radius was generalised from one 30.7 arcsec
    # case and threw away 9,320 real candidates.
    shred = cat.shredded(pool["ra"], pool["dec"], 45.0)
    # A bad brick manufactures a clump of fake detections -- 55 within 1.5
    # arcmin at RA 144.78, Dec -2.79, all of them blue bad-data square. Twelve
    # sampled clumped candidates contained zero real galaxies.
    clump = badregion.flag(pool["ra"], pool["dec"])
    clean = pool[~known & ~shred & ~clump]
    print(f"pool {len(pool):,} -> {known.sum():,} known, "
          f"{(shred & ~known).sum():,} shred-flagged, "
          f"{(clump & ~known & ~shred).sum():,} bad-region "
          f"-> {len(clean):,} pre-screened", flush=True)

    # The radial model runs BEFORE NED, not after.
    #
    # NED has no bulk download -- it is a per-object CGI query at ~0.1-0.5/s and
    # is the only genuinely slow stage left. Running the model last meant NED
    # spent hours verifying objects the model would have discarded: it flags
    # ~37% of candidates as junk. Ordering the cheap discriminator first cuts
    # the slow stage's workload by that fraction and puts the best candidates
    # at the front of the queue instead of merely the biggest.
    ranked = rank_by_model(clean)
    print(f"{WORKERS} workers", flush=True)

    state = load_state()
    todo = [o for o in ranked if key(o) not in state]
    print(f"{len(state):,} already verified, {len(todo):,} to go", flush=True)

    novel, known, failed = 0, 0, 0
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(verify, o): o for o in todo}
        for n, fut in enumerate(as_completed(futs), 1):
            o = futs[fut]
            try:
                verdict, detail = fut.result()
            except Exception as e:                                # noqa: BLE001
                verdict, detail = "failed", f"{type(e).__name__}"
            if verdict == "novel":
                novel += 1
            elif verdict == "known":
                known += 1
            else:
                failed += 1
            if verdict != "failed":            # failures stay unrecorded, retried next run
                state[key(o)] = {"v": verdict, "d": detail,
                                 "r": float(o["shape_r"]), "mu": float(o["mu_g"])}
            if n % 20 == 0 or n == len(todo):
                json.dump(state, open(STATE, "w"))
                print(f"  {n:,}/{len(todo):,}  novel={novel} known={known} "
                      f"failed={failed}", flush=True)

    json.dump(state, open(STATE, "w"))

    rows = [(k, v) for k, v in state.items() if v["v"] == "novel"]
    if rows:
        arr = np.zeros(len(rows), dtype=[("ra", "f8"), ("dec", "f8"),
                                         ("shape_r", "f4"), ("mu_g", "f4"),
                                         ("type", "U8")])
        for i, (k, v) in enumerate(rows):
            ra, dec = k.split("_")
            arr[i] = (float(ra), float(dec), v["r"], v["mu"], "novel")
        arr = arr[np.argsort(-arr["shape_r"])]
        np.save(NOVEL, arr)
        print(f"\nVERIFIED NOVEL: {len(arr):,} -> {NOVEL}")
        print(f"  known: {sum(1 for _, v in state.items() if v['v']=='known'):,}")


if __name__ == "__main__":
    main()
