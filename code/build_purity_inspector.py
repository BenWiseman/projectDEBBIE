"""Build the purity inspector for the v1.2 stratified sample.

Why this exists as a second file rather than an edit to build_inspector.py:

  build_inspector.py builds inspect.html from a UNIFORM draw over
  catalogue_winnie/v1.1_main.npy. That file was written at 17:40 on 4 September;
  the sample the paper actually reserves for purity was drawn at 20:03 the same
  day, from v1.2 with the DESI-identified entries removed. So inspect.html
  serves the wrong catalogue: only 1,817 of its 2,000 objects are still in
  v1.2. Inspecting it would measure the purity of a catalogue we are not
  publishing. The old builder is left alone because the release cites it.

A third app, purity.html, already served this same v2 sample and weighted it
correctly; it is retired in favour of this one because it used a normal
approximation for the interval instead of Wilson at the effective sample size,
showed no per-cell progress, and had no builder in the repository, so it could
not be regenerated. Its localStorage key is deliberately not reused: see KEY.

Two things this must get right and the old app does not:

  1. THE SAMPLE. inspect_sample_v2.json, seed 20260905, 1,999 objects over 28
     cells of the (r_eff, mu) plane, allocated 70% proportional / 30% equal.
  2. THE ARITHMETIC. The draw is not uniform, so an unweighted yes/(yes+no) is
     biased towards the over-sampled large-r cells. Every object carries the
     weight of the catalogue entries it stands for: one object in the smallest
     cell stands for about 3 entries, one in r3-3.5 mu26-26.51 for about 106.
     Purity is the weight-reconstructed proportion, and its interval uses the
     effective sample size (1,775.6, not 1,999) because the weights are
     unequal. A naive binomial over the object count would quote an interval
     about 6% too narrow.

The browser's number is a progress display, not the result. Every verdict line
carries its weight and cell so purity_from_verdicts.py can recompute the figure
outside the browser, which is what the paper should quote.

    .venv/bin/python sweeps/build_purity_inspector.py
"""
from __future__ import annotations

import io
import json
import os

LEGACY = "/home/ben/legacy"
SAMPLE = os.path.join(LEGACY, "inspect_sample_v2.json")
OUT = os.path.join(LEGACY, "inspect_v2.html")
CATALOGUE_N = 161721          # v1.2 main, |b| >= 30, DESI-identified removed


def main():
    meta = json.load(io.open(SAMPLE, encoding="utf-8"))
    objs = meta["objects"]

    missing = [o for o in objs if not os.path.exists(os.path.join(LEGACY, o["img"]))]
    if missing:
        raise SystemExit(f"{len(missing)} of {len(objs)} cutouts are not on disk; "
                         f"first missing: {missing[0]['img']}")

    wsum = sum(o["w"] for o in objs)
    neff = wsum * wsum / sum(o["w"] ** 2 for o in objs)
    print(f"sample   {len(objs):,} objects, {len({o['cell'] for o in objs})} cells, "
          f"seed {meta['seed']}")
    print(f"source   {meta['source']}")
    print(f"weights  sum {wsum:,.0f} against catalogue {CATALOGUE_N:,} "
          f"({wsum / CATALOGUE_N:.4%}); n_eff {neff:,.1f}")
    print(f"images   {len(objs):,}/{len(objs):,} present")

    payload = json.dumps([
        {"ra": o["ra"], "dec": o["dec"], "r": o["r"], "mu": o["mu"],
         "w": o["w"], "cell": o["cell"], "img": o["img"]}
        for o in objs
    ], separators=(",", ":"))

    doc = TEMPLATE.replace("__OBJECTS__", payload) \
                  .replace("__NEFF__", f"{neff:.1f}") \
                  .replace("__WSUM__", f"{wsum:.0f}") \
                  .replace("__SEED__", str(meta["seed"])) \
                  .replace("__NOBJ__", str(len(objs)))
    io.open(OUT, "w", encoding="utf-8").write(doc)
    print(f"\nwrote {OUT}")
    print(f"  open with:  xdg-open {OUT}")


TEMPLATE = r"""<!doctype html><meta charset=utf-8><title>DEBBIE purity inspector</title>
<style>
:root{color-scheme:dark}
*{box-sizing:border-box}
body{background:#0b0b0e;color:#ddd;font:14px/1.5 system-ui,sans-serif;margin:0;
     display:flex;flex-direction:column;height:100vh;overflow:hidden}
header{padding:8px 14px;border-bottom:1px solid #222;display:flex;gap:16px;
       align-items:center;flex-wrap:wrap}
header b{color:#eee}
#bar{flex:1;height:6px;background:#1a1a20;border-radius:3px;overflow:hidden;min-width:110px}
#fill{height:100%;background:#4a9;width:0;transition:width .15s}
main{flex:1;display:flex;align-items:center;justify-content:center;position:relative;min-height:0}
#img{max-height:76vh;max-width:76vw;border-radius:6px;background:#000}
#meta{position:absolute;bottom:6px;left:0;right:0;text-align:center;color:#9aa;font-size:12px}
#purity{position:absolute;top:8px;right:12px;text-align:right;font-size:12px;color:#9aa;
        background:#0b0b0ecc;padding:6px 9px;border-radius:5px;line-height:1.45}
#purity b{color:#dde;font-size:15px}
#cells{position:absolute;top:8px;left:12px;font-size:11px;color:#778;
       background:#0b0b0ecc;padding:6px 9px;border-radius:5px;max-height:78vh;overflow:auto}
#cells div{white-space:nowrap}
#cells span{color:#5c8}
footer{padding:8px 14px;border-top:1px solid #222;color:#889;font-size:12px;
       display:flex;gap:18px;justify-content:center;flex-wrap:wrap}
kbd{background:#1c1c22;border:1px solid #333;border-radius:4px;padding:1px 6px;color:#cde}
#flash{position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);font-size:64px;
       font-weight:700;opacity:0;pointer-events:none;transition:opacity .12s}
#done{display:none;padding:32px;text-align:center;line-height:1.9;max-width:680px}
#done h2{margin:0 0 6px}
.note{color:#778;font-size:12px}
</style>
<header>
  <b>DEBBIE purity</b>
  <span id=pos></span>
  <div id=bar><div id=fill></div></div>
  <span id=tally></span>
  <span class=note>seed __SEED__ &middot; n_eff __NEFF__</span>
</header>
<main>
  <img id=img alt="">
  <div id=flash></div>
  <div id=cells></div>
  <div id=purity></div>
  <div id=meta></div>
  <div id=done></div>
</main>
<footer>
  <span><kbd>&larr;</kbd> not a galaxy</span>
  <span><kbd>&rarr;</kbd> galaxy</span>
  <span><kbd>&uarr;</kbd> ambiguous</span>
  <span><kbd>&darr;</kbd> second opinion</span>
  <span><kbd>z</kbd> undo</span>
  <span><kbd>s</kbd> save</span>
</footer>
<script>
const OBJ  = __OBJECTS__;
const NEFF = __NEFF__, WSUM = __WSUM__, NOBJ = __NOBJ__;
// Keyed to the sample's seed, not to a version word. A pre-existing purity.html
// used "debbie_purity_v2" with a DIFFERENT verdict-id format (raw float concat
// rather than 5dp), so sharing that key between the two apps would have made each
// blind to the other's verdicts while both wrote under it: the union would
// double-count every object judged twice. That app is retired; this key cannot
// collide with anything it left behind.
const KEY  = "debbie_purity_sample_20260905";
let verdicts = {};
try { verdicts = JSON.parse(localStorage.getItem(KEY) || "{}"); } catch (e) { verdicts = {}; }
let i = 0;
const $ = id => document.getElementById(id);
const img=$("img"), pos=$("pos"), fill=$("fill"), meta=$("meta"), tally=$("tally"),
      flash=$("flash"), done=$("done"), purity=$("purity"), cells=$("cells");

const id = o => o.ra.toFixed(5) + "_" + o.dec.toFixed(5);

function firstUndone(){
  for (let k = 0; k < OBJ.length; k++) if (!(id(OBJ[k]) in verdicts)) return k;
  return OBJ.length;
}

/* Wilson interval on the WEIGHTED proportion, using the effective sample size of the
   decided subset. The weights are unequal, so a binomial over the object count would
   be too narrow. */
function wilson(p, n){
  if (!(n > 0)) return [0, 0];
  const z = 1.959964, z2 = z * z;
  const d = 1 + z2 / n;
  const c = (p + z2 / (2 * n)) / d;
  const h = (z / d) * Math.sqrt(p * (1 - p) / n + z2 / (4 * n * n));
  return [Math.max(0, c - h), Math.min(1, c + h)];
}

function stats(){
  let wy = 0, wn = 0, wa = 0, ws = 0, sw = 0, sw2 = 0, cy = 0, cn = 0, ca = 0, cs = 0;
  for (const v of Object.values(verdicts)){
    if (v.v === "yes")         { wy += v.w; cy++; }
    else if (v.v === "no")     { wn += v.w; cn++; }
    else if (v.v === "amb")    { wa += v.w; ca++; }
    else if (v.v === "second") { ws += v.w; cs++; }
    if (v.v === "yes" || v.v === "no"){ sw += v.w; sw2 += v.w * v.w; }
  }
  const dw = wy + wn;
  const p  = dw > 0 ? wy / dw : 0;
  const ne = sw2 > 0 ? (sw * sw) / sw2 : 0;      // effective n of the DECIDED subset
  return {wy, wn, wa, ws, cy, cn, ca, cs, p, ne, dw,
          ci: wilson(p, ne), done: cy + cn + ca + cs};
}

function drawCells(){
  const tot = {}, got = {};
  for (const o of OBJ) tot[o.cell] = (tot[o.cell] || 0) + 1;
  for (const v of Object.values(verdicts)) got[v.cell] = (got[v.cell] || 0) + 1;
  cells.innerHTML = Object.keys(tot).sort().map(c =>
    "<div>" + c.padEnd(22, " ") + " <span>" + (got[c] || 0) + "</span>/" + tot[c] + "</div>"
  ).join("");
}

function show(){
  if (i >= OBJ.length){ finish(); return; }
  const o = OBJ[i], s = stats();
  img.src = o.img; img.style.display = ""; done.style.display = "none";
  pos.textContent = (i + 1) + " / " + OBJ.length;
  fill.style.width = (100 * s.done / OBJ.length) + "%";
  meta.textContent = o.ra.toFixed(4) + " " + (o.dec >= 0 ? "+" : "") + o.dec.toFixed(4)
    + "   r_e=" + o.r.toFixed(1) + "″   μ=" + o.mu.toFixed(1)
    + "   " + o.cell + "   stands for " + o.w.toFixed(0) + " entries";
  tally.innerHTML = "<b>" + s.cy + "</b> yes &nbsp;<b>" + s.cn + "</b> no &nbsp;<b>"
    + s.ca + "</b> ambig &nbsp;<b>" + s.cs + "</b> 2nd";
  purity.innerHTML = s.dw > 0
    ? "weighted purity<br><b>" + (100 * s.p).toFixed(1) + "%</b> ["
      + (100 * s.ci[0]).toFixed(1) + ", " + (100 * s.ci[1]).toFixed(1) + "]<br>"
      + "n_eff decided " + s.ne.toFixed(0)
    : "weighted purity<br><b>&mdash;</b><br>no decisions yet";
  drawCells();
}

function finish(){
  const s = stats();
  img.style.display = "none"; meta.textContent = ""; purity.innerHTML = ""; cells.innerHTML = "";
  done.style.display = "block";
  done.innerHTML = "<h2>Sample complete</h2>"
    + "<p><b>" + s.cy + "</b> galaxies, <b>" + s.cn + "</b> not, <b>" + s.ca
    + "</b> ambiguous, <b>" + s.cs + "</b> flagged for a second look</p>"
    + "<p style='font-size:22px'>weighted purity <b>" + (100 * s.p).toFixed(1) + "%</b> "
    + "[" + (100 * s.ci[0]).toFixed(1) + ", " + (100 * s.ci[1]).toFixed(1) + "]</p>"
    + "<p class=note>Weight-reconstructed over " + s.dw.toFixed(0) + " catalogue entries; "
    + "Wilson interval on n_eff = " + s.ne.toFixed(0) + " for the decided subset, "
    + "not on the " + (s.cy + s.cn) + " objects. Ambiguous entries are excluded from "
    + "both numerator and denominator, which assumes they split like the decided ones.</p>"
    + "<p>Press <kbd>s</kbd> to save the verdicts, then recompute the figure with "
    + "<code>purity_from_verdicts.py</code>. The number above is a progress display; "
    + "the paper should quote the script.</p>";
}

function mark(v){
  if (i >= OBJ.length) return;
  const o = OBJ[i];
  verdicts[id(o)] = {v: v, ra: o.ra, dec: o.dec, r: o.r, mu: o.mu, w: o.w,
                     cell: o.cell, t: new Date().toISOString()};
  try { localStorage.setItem(KEY, JSON.stringify(verdicts)); } catch (e) {}
  flash.textContent = {yes:"✓", no:"✕", amb:"?", second:"⚠"}[v];
  flash.style.color = {yes:"#5c8", no:"#c55", amb:"#dc6", second:"#59c"}[v];
  flash.style.opacity = 1; setTimeout(() => flash.style.opacity = 0, 120);
  i++; show();
}

function undo(){
  if (i <= 0) return;
  i--; delete verdicts[id(OBJ[i])];
  try { localStorage.setItem(KEY, JSON.stringify(verdicts)); } catch (e) {}
  show();
}

function save(){
  const lines = Object.values(verdicts).map(v => JSON.stringify(v)).join("\n");
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([lines + "\n"], {type: "application/x-ndjson"}));
  a.download = "purity_verdicts_v2.jsonl";     // distinct from the shortlist's file
  a.click();
}

addEventListener("keydown", e => {
  const k = e.key;
  if (k === "ArrowLeft")       { mark("no");     e.preventDefault(); }
  else if (k === "ArrowRight") { mark("yes");    e.preventDefault(); }
  else if (k === "ArrowUp")    { mark("amb");    e.preventDefault(); }
  else if (k === "ArrowDown")  { mark("second"); e.preventDefault(); }
  else if (k === "z")          { undo(); }
  else if (k === "s")          { save(); }
});

function preload(){
  for (let k = i + 1; k < Math.min(i + 8, OBJ.length); k++){
    const im = new Image(); im.src = OBJ[k].img;
  }
}
setInterval(preload, 400);
i = firstUndone(); show();
</script>
"""


if __name__ == "__main__":
    main()
