"""Weighted catalogue purity from the inspector's verdict file.

The browser shows a running number so the inspection is not flown blind. This is
the number the paper quotes, because it can be re-run and diffed.

    .venv/bin/python sweeps/purity_from_verdicts.py ~/Downloads/purity_verdicts_v2.jsonl
    .venv/bin/python sweeps/purity_from_verdicts.py --selftest

The arithmetic, and why it is not yes/(yes+no):

  The sample is stratified 70% proportional / 30% equal over 28 cells of the
  (r_eff, mu) plane, so one object in the smallest cell stands for about 3
  catalogue entries and one in r3-3.5 mu26-26.51 stands for about 106. An
  unweighted proportion is therefore a proportion of the SAMPLE, not of the
  catalogue, and it is biased towards the deliberately over-sampled large-r
  cells. Each verdict carries its weight; purity is the weight-reconstructed
  proportion and the interval uses the effective sample size of the decided
  subset, n_eff = (sum w)^2 / sum(w^2).

Ambiguous verdicts are excluded from numerator and denominator alike, which
assumes they would split like the decided ones. That assumption is reported as a
bound: the printout gives purity with every ambiguous entry counted as a galaxy
and as a non-galaxy, and the truth is inside that range whatever they are.
"""
from __future__ import annotations

import collections
import io
import json
import math
import sys

Z = 1.959964          # two-sided 95%
CATALOGUE_N = 161721  # v1.2 main, |b| >= 30


def wilson(k_over_n, n):
    """Wilson score interval at an effective (possibly non-integer) sample size."""
    if n <= 0:
        return (0.0, 0.0)
    p, z2 = k_over_n, Z * Z
    d = 1.0 + z2 / n
    c = (p + z2 / (2 * n)) / d
    h = (Z / d) * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n))
    return (max(0.0, c - h), min(1.0, c + h))


def weighted_purity(rows):
    """Weight-reconstructed purity over the decided rows, with a Wilson interval on n_eff."""
    wy = sum(r["w"] for r in rows if r["v"] == "yes")
    wn = sum(r["w"] for r in rows if r["v"] == "no")
    dec = [r for r in rows if r["v"] in ("yes", "no")]
    sw = sum(r["w"] for r in dec)
    sw2 = sum(r["w"] ** 2 for r in dec)
    neff = (sw * sw / sw2) if sw2 > 0 else 0.0
    p = (wy / (wy + wn)) if (wy + wn) > 0 else 0.0
    return p, wilson(p, neff), neff, wy, wn


def report(path):
    rows = [json.loads(l) for l in io.open(path, encoding="utf-8") if l.strip()]
    for r in rows:
        if "w" not in r or "v" not in r:
            raise SystemExit(f"verdict line without a weight or verdict: {r}. "
                             f"This file did not come from build_purity_inspector.py.")

    by = collections.Counter(r["v"] for r in rows)
    p, ci, neff, wy, wn = weighted_purity(rows)

    print(f"verdicts        {len(rows):,}")
    for k in ("yes", "no", "amb", "second"):
        print(f"  {k:<12s}  {by.get(k, 0):,}")
    print(f"cells covered   {len({r['cell'] for r in rows})} of 28")
    print()
    print(f"WEIGHTED PURITY {100 * p:.2f}%  95% CI [{100 * ci[0]:.2f}, {100 * ci[1]:.2f}]")
    print(f"  reconstructed over {wy + wn:,.0f} catalogue entries of {CATALOGUE_N:,}")
    print(f"  n_eff (decided)    {neff:,.1f}")

    unweighted = by.get("yes", 0) / max(1, by.get("yes", 0) + by.get("no", 0))
    print(f"  unweighted, for comparison only: {100 * unweighted:.2f}% "
          f"(do NOT quote: a proportion of the sample, not of the catalogue)")

    # EVERY undecided verdict has to be bounded, not just the ambiguous ones.
    # "second opinion" is also undecided, and on the v1.2 sample it is the larger
    # of the two; bounding only "amb" understated the range badly.
    und = [r for r in rows if r["v"] not in ("yes", "no")]
    if und:
        decided = [r for r in rows if r["v"] in ("yes", "no")]
        lo, _, _, _, _ = weighted_purity(decided + [dict(r, v="no") for r in und])
        hi, _, _, _, _ = weighted_purity(decided + [dict(r, v="yes") for r in und])
        kinds = collections.Counter(r["v"] for r in und)
        print()
        print(f"  {len(und):,} undecided are excluded above "
              f"({', '.join(f'{v} {k}' for k, v in sorted(kinds.items()))}).")
        print(f"  Counting them all as non-galaxies gives {100 * lo:.2f}%, "
              f"all as galaxies {100 * hi:.2f}%.")
        print(f"  Purity is inside [{100 * lo:.2f}, {100 * hi:.2f}] whatever they are,")
        print(f"  and that range is {100 * (hi - lo):.1f} points wide against a "
              f"sampling interval of {100 * (ci[1] - ci[0]):.1f}.")

    if len(rows) < 1999:
        print()
        print(f"  INCOMPLETE: {1999 - len(rows):,} of 1,999 objects not yet judged. "
              f"Strata are unevenly filled, so this is not yet the sample's answer.")
    return 0


def selftest():
    """Plant cases with known answers, including ones that must FAIL if the weighting is dropped."""
    ok = True

    def check(name, got, want, tol=1e-6):
        nonlocal ok
        good = abs(got - want) < tol
        ok = ok and good
        print(f"  {name:<58s} {'ok' if good else 'FAIL got %.6f want %.6f' % (got, want)}")

    # Uniform weights: weighted must equal unweighted.
    rows = [{"v": "yes", "w": 1.0}] * 8 + [{"v": "no", "w": 1.0}] * 2
    p, ci, neff, _, _ = weighted_purity(rows)
    check("uniform weights reproduce the unweighted proportion", p, 0.8)
    check("uniform weights give n_eff = n", neff, 10.0)

    # The case the old app got wrong: many heavy 'no', few light 'yes'.
    # Unweighted says 50%; the catalogue says 9.09%.
    rows = [{"v": "yes", "w": 1.0}] * 10 + [{"v": "no", "w": 10.0}] * 10
    p, _, neff, _, _ = weighted_purity(rows)
    check("unequal weights do NOT give the unweighted 50%", p, 10 / 110.0)
    check("  and n_eff is below the object count", 1.0 if neff < 20 else 0.0, 1.0)

    # n_eff formula on a known asymmetric case.
    rows = [{"v": "yes", "w": 3.0}, {"v": "no", "w": 1.0}]
    _, _, neff, _, _ = weighted_purity(rows)
    check("n_eff = (sum w)^2 / sum w^2 = 16/10", neff, 1.6)

    # Wilson at an effective n must be wider than at the object count when n_eff < n.
    _, ci_eff = wilson(0.9, 1775.6), None
    w_eff = wilson(0.9, 1775.6)
    w_obj = wilson(0.9, 1999.0)
    wider = (w_eff[1] - w_eff[0]) > (w_obj[1] - w_obj[0])
    check("Wilson on n_eff is wider than on the object count", 1.0 if wider else 0.0, 1.0)

    # A planted failure: if someone reintroduces an unweighted count, this must catch it.
    def unweighted(rows):
        y = sum(1 for r in rows if r["v"] == "yes")
        n = sum(1 for r in rows if r["v"] == "no")
        return y / (y + n)
    rows = [{"v": "yes", "w": 1.0}] * 10 + [{"v": "no", "w": 10.0}] * 10
    caught = abs(unweighted(rows) - weighted_purity(rows)[0]) > 0.1
    check("planted: an unweighted count is detectably different here",
          1.0 if caught else 0.0, 1.0)

    print("\n  ALL PASS" if ok else "\n  FAILURES ABOVE")
    return 0 if ok else 1


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(selftest())
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    sys.exit(report(sys.argv[1]))
