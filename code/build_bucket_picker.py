"""Build the bucket exemplar picker: one representative image per (r_eff, mu) cell.

The figure this feeds is a grid of real catalogue objects laid out on the two axes
the selection is defined on, so a reader can see what a 3.2 arcsec object at
mu 26.3 actually looks like against a 12 arcsec object at mu 24.5.

## Two groups per bucket, and why they are not interchangeable

Each bucket offers CONFIRMED objects (ones a human has judged a galaxy, from the
385-object shortlist pass and from the purity sample) and TYPICAL objects (a
seeded random draw from every catalogue entry in that bucket). They support
different captions and must not be mixed silently: a confirmed panel says "this
is what a real one looks like here", a typical panel says "this is what the
catalogue contains here". The export records which group each pick came from.

The shortlist cut at r_eff >= 6 arcsec, so it only populates the 12 large-radius
buckets. The 16 small-radius buckets fill in from the purity verdicts, which span
all 28 cells. Re-run this after exporting purity_verdicts_v2.jsonl.

## The honesty problem, and how this handles it

Hand-picking the best-looking object in each cell produces a figure that
misrepresents the catalogue, and it is the easiest mistake in the world to make
because the nice ones are the ones you want to show. Equally, a blind random
draw will sometimes land on a diffraction spike and waste a cell.

So: each cell is served as a SEEDED RANDOM draw, always in the same order, and
the export records how deep into that order the chosen tile sat. A caption can
then say "each panel is the first acceptable object from a seeded random draw,
median depth N" and a referee can check it. Picking the 1st is free; picking the
19th is visible.

Cells are the same 28 used by the purity sample, so the figure, the selection
function and the purity table all cut the plane the same way.

    .venv/bin/python sweeps/build_bucket_picker.py [n_per_bucket]
"""
from __future__ import annotations

import io
import json
import os
import sys

import numpy as np

LEGACY = "/home/ben/legacy"
ARCH = os.path.join(LEGACY, "candidate_images")
CAT = "/home/ben/star/catalogue_debbie/v1.2_main.npy"
SHORTLIST = "/home/ben/star/catalogue_debbie/shortlist_verdicts_2026-09-30.csv"
PURITY = os.path.expanduser("~/Downloads/purity_verdicts_v2.jsonl")
OUT = os.path.join(LEGACY, "buckets.html")
SEED = 20261001
N_PER = int(sys.argv[1]) if len(sys.argv) > 1 else 24

R_EDGES = [3, 3.5, 4, 5, 6, 8, 10, 15.1]
M_EDGES = [24, 25, 25.5, 26, 26.51]


def img_rel(ra, dec):
    return f"candidate_images/ra{int(ra):03d}/{ra:.5f}{dec:+.5f}_75.jpg"


def main():
    cat = np.load(CAT, allow_pickle=True)
    ra = np.asarray(cat["ra"], float)
    dec = np.asarray(cat["dec"], float)
    r = np.asarray(cat["shape_r"], float)
    mu = np.asarray(cat["mu_g"], float)
    print(f"catalogue {len(ra):,}")

    # every human verdict available, from both passes
    import csv
    verdict, source = {}, {}
    if os.path.exists(SHORTLIST):
        n = 0
        for row in csv.DictReader(io.open(SHORTLIST, encoding="utf-8")):
            k = (round(float(row["ra_deg"]), 5), round(float(row["dec_deg"]), 5))
            verdict[k] = row["verdict"]; source[k] = "shortlist"; n += 1
        print(f"shortlist verdicts  {n:,}")
    if os.path.exists(PURITY):
        raw = io.open(PURITY, encoding="utf-8").read()
        # an earlier exporter wrote literal backslash-n; tolerate both forms
        lines = [l for l in raw.replace(chr(92) + "n", chr(10)).split(chr(10)) if l.strip()]
        for line in lines:
            v = json.loads(line)
            k = (round(float(v["ra"]), 5), round(float(v["dec"]), 5))
            verdict.setdefault(k, v["v"]); source.setdefault(k, "purity")
        print(f"purity verdicts     {len(lines):,}")
    else:
        print(f"purity verdicts     none yet; export and re-run to fill the small-r buckets")
    print(f"confirmed total     {sum(1 for v in verdict.values() if v == 'yes'):,}")

    rng = np.random.default_rng(SEED)
    buckets, missing = [], 0
    for i in range(len(R_EDGES) - 1):
        for j in range(len(M_EDGES) - 1):
            lo_r, hi_r, lo_m, hi_m = R_EDGES[i], R_EDGES[i + 1], M_EDGES[j], M_EDGES[j + 1]
            sel = np.where((r >= lo_r) & (r < hi_r) & (mu >= lo_m) & (mu < hi_m))[0]
            if not len(sel):
                continue
            def tile(k, rank, group):
                key = (round(float(ra[k]), 5), round(float(dec[k]), 5))
                return {"ra": key[0], "dec": key[1],
                        "r": round(float(r[k]), 2), "mu": round(float(mu[k]), 2),
                        "img": img_rel(ra[k], dec[k]), "rank": rank, "group": group,
                        "v": verdict.get(key, ""), "src": source.get(key, "")}

            def have(k):
                return os.path.exists(os.path.join(LEGACY, img_rel(ra[k], dec[k])))

            # CONFIRMED and REJECTED: what a human kept and what they threw out.
            # Shown side by side, these are the figure that illustrates the screen.
            conf, rej = [], []
            for k in sel:
                k = int(k)
                v = verdict.get((round(float(ra[k]), 5), round(float(dec[k]), 5)))
                if v not in ("yes", "no"):
                    continue
                if not have(k):
                    missing += 1
                    continue
                if v == "yes":
                    conf.append(tile(k, len(conf) + 1, "confirmed"))
                else:
                    rej.append(tile(k, len(rej) + 1, "rejected"))

            # TYPICAL: a seeded random draw over the whole cell, order recorded
            order = rng.permutation(len(sel))
            tiles = []
            for pos in order:
                k = int(sel[pos])
                if not have(k):
                    missing += 1
                    continue
                tiles.append(tile(k, len(tiles) + 1, "typical"))
                if len(tiles) >= N_PER:
                    break
            buckets.append({"r_lo": lo_r, "r_hi": hi_r, "mu_lo": lo_m, "mu_hi": hi_m,
                            "n_cat": int(len(sel)), "confirmed": conf, "rejected": rej,
                            "typical": tiles, "n_conf": len(conf), "n_rej": len(rej),
                            "label": f"r {lo_r}-{hi_r}″   μ {lo_m}-{hi_m}"})
    n_conf = sum(b["n_conf"] for b in buckets)
    n_rej = sum(b["n_rej"] for b in buckets)
    n_typ = sum(len(b["typical"]) for b in buckets)
    print(f"rejected tiles      {n_rej}")
    empty = [b["label"] for b in buckets if not b["n_conf"]]
    print(f"\nbuckets {len(buckets)}: {n_conf} confirmed tiles, {n_typ} typical tiles, "
          f"{missing} cutouts missing")
    if empty:
        print(f"buckets with no confirmed object yet: {len(empty)} of {len(buckets)}")

    doc = (TEMPLATE
           .replace("__BUCKETS__", json.dumps(buckets, separators=(",", ":")))
           .replace("__SEED__", str(SEED))
           .replace("__NPER__", str(N_PER)))
    io.open(OUT, "w", encoding="utf-8").write(doc)
    print(f"\nwrote {OUT}")
    print(f"  open with:  xdg-open {OUT}")


TEMPLATE = r"""<!doctype html><meta charset=utf-8><title>DEBBIE bucket exemplars</title>
<style>
:root{color-scheme:dark}
*{box-sizing:border-box}
body{background:#0b0b0e;color:#ddd;font:14px/1.5 system-ui,sans-serif;margin:0;
     display:flex;flex-direction:column;height:100vh;overflow:hidden}
header{padding:8px 14px;border-bottom:1px solid #222;display:flex;gap:16px;align-items:center;flex-wrap:wrap}
header b{color:#eee}
.note{color:#778;font-size:12px}
main{flex:1;display:flex;min-height:0}
#side{width:210px;border-right:1px solid #222;overflow:auto;padding:6px 0;flex:none}
#side div{padding:4px 10px;cursor:pointer;font-size:12px;white-space:nowrap}
#side div:hover{background:#15151a}
#side div.on{background:#1d2a24;color:#8fd}
#side div.cur{background:#243; color:#cfe}
#side span{color:#5c8}
#grid{flex:1;overflow:auto;padding:10px;display:grid;
      grid-template-columns:repeat(auto-fill,minmax(140px,1fr));gap:8px;align-content:start}
.hdr{grid-column:1/-1;color:#8ab;font-size:11px;letter-spacing:.04em;
     padding:10px 2px 2px;border-bottom:1px solid #223;margin-bottom:4px}
figure{margin:0;cursor:pointer;border:2px solid transparent;border-radius:6px;padding:3px;position:relative}
figure:hover{border-color:#456}
figure.pick{border-color:#5c8;background:#13201a}
figure img{width:100%;display:block;border-radius:4px;background:#000}
figcaption{font-size:10px;color:#889;text-align:center;padding-top:2px}
.rank{position:absolute;top:5px;left:5px;background:#000a;color:#9ab;font-size:10px;
      padding:0 4px;border-radius:3px}
.vd{position:absolute;top:5px;right:5px;font-size:11px;padding:0 4px;border-radius:3px;background:#000a}
footer{padding:7px 14px;border-top:1px solid #222;color:#889;font-size:12px;
       display:flex;gap:18px;justify-content:center;flex-wrap:wrap}
kbd{background:#1c1c22;border:1px solid #333;border-radius:4px;padding:1px 6px;color:#cde}
</style>
<header>
  <b>DEBBIE bucket exemplars</b>
  <span id=pos></span>
  <span id=tally class=note></span>
  <span class=note>seed __SEED__ &middot; __NPER__ per bucket, in seeded random order</span>
</header>
<main>
  <div id=side></div>
  <div id=grid></div>
</main>
<footer>
  <span><kbd>&larr;</kbd><kbd>&rarr;</kbd> bucket</span>
  <span>click a tile to pick it for that bucket</span>
  <span><kbd>x</kbd> clear this bucket</span>
  <span><kbd>s</kbd> save picks</span>
</footer>
<script>
const B = __BUCKETS__;
const KEY = "debbie_bucket_picks_v1";
let picks = {};
try { picks = JSON.parse(localStorage.getItem(KEY) || "{}"); } catch(e) { picks = {}; }
let bi = 0;
const $ = id => document.getElementById(id);

function bkey(b){ return b.r_lo+"_"+b.mu_lo; }

function drawSide(){
  $("side").innerHTML = B.map((b,i) =>
    "<div class='"+(i===bi?"cur ":"")+(picks[bkey(b)]?"on":"")+"' data-i='"+i+"'>"
    + b.label + " <span>" + (picks[bkey(b)] ? "✓" : "") + "</span></div>").join("");
  [...$("side").children].forEach(d =>
    d.onclick = () => { bi = +d.dataset.i; draw(); });
}

function draw(){
  const b = B[bi], k = bkey(b), chosen = picks[k];
  $("pos").textContent = (bi+1) + " / " + B.length + "   " + b.label
    + "   (" + b.n_cat.toLocaleString() + " in catalogue)";
  $("tally").textContent = Object.keys(picks).length + " of " + B.length + " buckets picked";
  const tileHTML = t => {
    const on = chosen && chosen.ra === t.ra && chosen.dec === t.dec;
    const vd = t.v === "yes" ? "<span class=vd style='color:#5c8'>✓</span>"
             : t.v === "no" ? "<span class=vd style='color:#c55'>✕</span>"
             : t.v ? "<span class=vd style='color:#dc6'>?</span>" : "";
    return "<figure class='"+(on?"pick":"")+"'>"
      + "<span class=rank>#"+t.rank+"</span>" + vd
      + "<img loading=lazy src='"+t.img+"'>"
      + "<figcaption>r "+t.r.toFixed(1)+"″  μ "+t.mu.toFixed(1)+"</figcaption></figure>";
  };
  const all = b.confirmed.concat(b.rejected, b.typical);
  const head = (txt,n) => "<div class=hdr>"+txt+" ("+n+")</div>";
  $("grid").innerHTML =
      head("CONFIRMED by eye — kept by the screen", b.confirmed.length)
    + b.confirmed.map(tileHTML).join("")
    + head("REJECTED by eye — thrown out by the screen, the other half of the comparison",
           b.rejected.length)
    + b.rejected.map(tileHTML).join("")
    + head("TYPICAL — seeded random draw over all "+b.n_cat.toLocaleString()
           +" in this cell, screened and unscreened alike", b.typical.length)
    + b.typical.map(tileHTML).join("");
  [...$("grid").querySelectorAll("figure")].forEach((f,idx) => f.onclick = () => {
    const t = all[idx];
    if (chosen && chosen.ra === t.ra && chosen.dec === t.dec) delete picks[k];
    else picks[k] = {ra:t.ra, dec:t.dec, r:t.r, mu:t.mu, img:t.img, rank:t.rank,
                     group:t.group, v:t.v, src:t.src, label:b.label,
                     shown_confirmed:b.confirmed.length, shown_typical:b.typical.length,
                     r_lo:b.r_lo, r_hi:b.r_hi, mu_lo:b.mu_lo, mu_hi:b.mu_hi};
    try { localStorage.setItem(KEY, JSON.stringify(picks)); } catch(e) {}
    drawSide(); draw();
  });
  drawSide();
}

function save(){
  const rows = B.filter(b => picks[bkey(b)]).map(b => picks[bkey(b)]);
  const typ = rows.filter(r => r.group === "typical").map(r => r.rank).sort((a,c)=>a-c);
  const med = typ.length ? typ[Math.floor(typ.length/2)] : 0;
  const blob = {seed: __SEED__, n_per_bucket: __NPER__, buckets: B.length,
                picked: rows.length,
                from_confirmed: rows.filter(r => r.group === "confirmed").length,
                from_typical: typ.length,
                median_pick_depth_typical: med, picks: rows};
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([JSON.stringify(blob, null, 1)], {type:"application/json"}));
  a.download = "bucket_picks.json"; a.click();
}

addEventListener("keydown", e => {
  if (e.key === "ArrowRight"){ bi = Math.min(B.length-1, bi+1); draw(); e.preventDefault(); }
  else if (e.key === "ArrowLeft"){ bi = Math.max(0, bi-1); draw(); e.preventDefault(); }
  else if (e.key === "x"){ delete picks[bkey(B[bi])];
    try { localStorage.setItem(KEY, JSON.stringify(picks)); } catch(e) {} draw(); }
  else if (e.key === "s"){ save(); }
});
draw();
</script>
"""


if __name__ == "__main__":
    main()
