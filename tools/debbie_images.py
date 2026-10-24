#!/usr/bin/env python3
"""Find DEBBIE entries and download their Legacy Surveys DR10 images.

For each entry chosen, three cutouts centred on the catalogued position: the
DR10 image, the Tractor model, and the residual (image minus model), as JPEG
or as FITS in g, r and z. A gallery page, index.html, shows them side by side
with each entry's size, surface brightness and score, and links to the Legacy
Surveys viewer.

Choosing entries (combine freely; applied in this order)
  --name "DEBBIE 3" ...       by designation, or --number 3 5 ...
  --cone RA DEC ARCMIN        within a radius of a position (degrees, arcmin)
  --part extension            one part: original, recovered or extension
  --where "shape_r>=8" ...    any column, with < <= > >= == !=; rows with no
                              value (nan) never pass <, <=, >, >= or ==
  --sample N                  a random N of what is left (--seed to repeat)
  --stratify COL:K ...        spread the sample evenly over K quantile bins of
                              each COL, or COL=e0,e1,...,en for fixed edges
  --sort COL [--desc]         order by a column, largest first with --desc
  --limit N                   keep the first N

Examples
  python3 tools/debbie_images.py --sample 60 --stratify deepshadows_mean:4 shape_r:3
  python3 tools/debbie_images.py --where "shape_r>=12" --sort deepshadows_mean --desc --limit 40
  python3 tools/debbie_images.py --part extension --sample 30 --context
  python3 tools/debbie_images.py --cone 185.45 -8.85 20 --format fits
  python3 tools/debbie_images.py --name "DEBBIE 3" "DEBBIE 5"
  python3 tools/debbie_images.py --sample 500 --stratify mu_eff_g:5 --dry-run

Output (--out, default debbie_images/):
  DEBBIE_3/image.jpg, model.jpg, resid.jpg (and context.jpg with --context)
  selection.csv   the chosen rows, with "cell" and "weight" first; weight is
                  how many selected entries each one stands for, so weighted
                  shares over a stratified sample estimate the whole selection
  index.html      the gallery
  failures.csv    any cutout that could not be fetched
Files already present are kept, so a stopped run carries on where it left off.

The field is max(--fov-min, --fov-reff x SHAPE_R) at the DR10 pixel scale,
0.262"; --context adds an image four times wider. The cutout service is
shared: the default is 4 requests at a time (at most 8) and no more than
1,000 entries a run unless --max-objects says otherwise.

Standard library only; Python 3.8 or later.
"""
import argparse
import bisect
import csv
import gzip
import html
import math
import os
import random
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
CATALOGUE = os.path.join(HERE, "..", "catalogue", "DEBBIE_v1.5.csv.gz")
SERVICE = "https://www.legacysurvey.org/viewer"
PIXSCALE = 0.262
LAYERS = {"image": "ls-dr10", "model": "ls-dr10-model", "resid": "ls-dr10-resid"}
MAGIC = {"jpg": b"\xff\xd8\xff", "fits": b"SIMPLE  ="}
WAITS = (0, 5, 15, 30)          # seconds before each try; the service sometimes answers "busy"
USER_AGENT = "debbie-images/1.0 (DEBBIE catalogue cutouts)"
SHOWN = ("part", "type", "shape_r", "mu_eff_g", "deepshadows_mean", "centroid_flag")
OPS = {"<": lambda a, b: a < b, "<=": lambda a, b: a <= b, ">": lambda a, b: a > b,
       ">=": lambda a, b: a >= b, "==": lambda a, b: a == b, "!=": lambda a, b: a != b}
WHERE = re.compile(r"^\s*(\w+)\s*(<=|>=|==|!=|<|>)\s*(.+?)\s*$")


class Catalogue:
    """The table as raw lines; columns are parsed when first asked for."""

    def __init__(self, path):
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", newline="") as f:
            self.header = f.readline().rstrip("\r\n").split(",")
            self.lines = [line.rstrip("\r\n") for line in f]
        if any('"' in line for line in self.lines):
            raise SystemExit(f"{path} has quoted fields; this reader splits on commas")
        self.index = {c: i for i, c in enumerate(self.header)}
        self.cache = {}

    def __len__(self):
        return len(self.lines)

    def column(self, name):
        if name not in self.index:
            raise SystemExit(f"no column {name!r}. Columns: {', '.join(self.header)}")
        if name not in self.cache:
            i = self.index[name]
            raw = [line.split(",", i + 1)[i] for line in self.lines]
            try:
                self.cache[name] = [float(v) if v else math.nan for v in raw]
            except ValueError:
                self.cache[name] = raw
        return self.cache[name]

    def numeric(self, name):
        col = self.column(name)
        if col and not isinstance(col[0], float):
            raise SystemExit(f"{name} is not a numeric column")
        return col

    def row(self, i):
        return self.lines[i].split(",")


def number_of(name):
    m = re.fullmatch(r"\s*(?:DEBBIE[\s_]*)?0*(\d+)\s*", str(name), re.I)
    if not m:
        raise SystemExit(f"not a DEBBIE designation: {name!r}")
    return int(m.group(1))


def by_number(cat, idx, numbers):
    want = set(numbers)
    num = cat.numeric("number")
    out = [i for i in idx if int(num[i]) in want]
    found = {int(num[i]) for i in out}
    lost = sorted(want - found)
    if lost:
        print(f"not in the catalogue (retired or never assigned): "
              f"{', '.join(f'DEBBIE {n}' for n in lost[:20])}{' ...' if len(lost) > 20 else ''}")
    return out


def separation(ra1, dec1, ra2, dec2):
    """Angle between two positions, all in degrees."""
    p1, p2 = math.radians(dec1), math.radians(dec2)
    h = (math.sin((p2 - p1) / 2) ** 2
         + math.cos(p1) * math.cos(p2) * math.sin(math.radians(ra2 - ra1) / 2) ** 2)
    return math.degrees(2 * math.asin(min(1.0, math.sqrt(h))))


def cone(cat, idx, ra0, dec0, arcmin):
    ra, dec = cat.numeric("ra"), cat.numeric("dec")
    r = arcmin / 60
    return [i for i in idx if abs(dec[i] - dec0) <= r and separation(ra0, dec0, ra[i], dec[i]) <= r]


def where(cat, idx, expr):
    m = WHERE.match(expr)
    if not m:
        raise SystemExit(f"cannot read --where {expr!r}; write it like shape_r>=8 or part==extension")
    name, op, value = m.groups()
    col = cat.column(name)
    if col and isinstance(col[0], float):
        try:
            value = float(value)
        except ValueError:
            raise SystemExit(f"{name} is numeric; {value!r} is not a number")
    elif op not in ("==", "!="):
        raise SystemExit(f"{name} is text; use == or !=")
    test = OPS[op]
    return [i for i in idx if test(col[i], value)]


def bins(cat, idx, spec):
    """(name, edges, kept rows) for one --stratify spec, COL:K or COL=e0,...,en."""
    m = re.fullmatch(r"(\w+)(?::(\d+)|=([-\d.,eE+]+))", spec)
    if not m:
        raise SystemExit(f"cannot read --stratify {spec!r}; write it like shape_r:3 or shape_r=3,5,8,20")
    name, k, edges = m.groups()
    col = cat.numeric(name)
    kept = [i for i in idx if not math.isnan(col[i])]
    if len(kept) < len(idx):
        print(f"{len(idx) - len(kept):,} rows with no {name} left out of the stratified sample")
    if edges:
        e = sorted(float(x) for x in edges.split(","))
        inside = [i for i in kept if e[0] <= col[i] <= e[-1]]
        if len(inside) < len(kept):
            print(f"{len(kept) - len(inside):,} rows with {name} outside {e[0]}-{e[-1]} left out")
        return name, e, inside
    k = int(k)
    vals = sorted(col[i] for i in kept)
    if not vals:
        return name, [0.0, 0.0], kept
    e = sorted(set([vals[0]] + [vals[min(len(vals) - 1, (j * len(vals)) // k)] for j in range(1, k)]
                   + [vals[-1]]))
    return name, e, kept


def cell_of(value, edges):
    """Bin number for value; the top edge is inside the last bin."""
    return min(len(edges) - 2, max(0, bisect.bisect_right(edges, value) - 1))


def stratified(cat, idx, specs, n, rng):
    """A sample of n spread as evenly as the cells allow, with its cells and weights."""
    cols = []
    for spec in specs:
        name, edges, idx = bins(cat, idx, spec)
        cols.append((name, edges))
    cells = {}
    for i in idx:
        key = tuple(cell_of(cat.numeric(name)[i], edges) for name, edges in cols)
        cells.setdefault(key, []).append(i)
    share = {key: 0 for key in cells}
    left = n
    while left > 0:
        room = [key for key in cells if share[key] < len(cells[key])]
        if not room:
            break
        rng.shuffle(room)
        each = max(1, left // len(room))
        for key in room:
            take = min(each, len(cells[key]) - share[key], left)
            share[key] += take
            left -= take
            if left == 0:
                break
    out, cell, weight = [], {}, {}
    for key in sorted(cells):
        if not share[key]:
            continue
        label = " ".join(f"{edges[b]:.4g}<={name}{'<=' if b == len(edges) - 2 else '<'}{edges[b + 1]:.4g}"
                         for (name, edges), b in zip(cols, key))
        for i in rng.sample(cells[key], share[key]):
            out.append(i)
            cell[i] = label
            weight[i] = len(cells[key]) / share[key]
    print(f"stratified: {len(out):,} from {sum(1 for k in share if share[k])} of {len(cells)} cells")
    return out, cell, weight


def order(cat, idx, key, desc=False):
    """Sorted by a numeric column; rows with no value go last."""
    col = cat.numeric(key)
    good = sorted((i for i in idx if not math.isnan(col[i])), key=lambda i: col[i], reverse=desc)
    return good + [i for i in idx if math.isnan(col[i])]


def field(shape_r, a):
    """(size in pixels, pixel scale) for the main cutouts."""
    fov = max(a.fov_min, a.fov_reff * shape_r)
    return max(64, min(3000, math.ceil(fov / PIXSCALE))), PIXSCALE


def cutout_url(layer, ra, dec, size, pixscale, fmt):
    q = {"ra": f"{ra:.7f}", "dec": f"{dec:.7f}", "layer": layer,
         "pixscale": f"{pixscale:.4f}", "size": size}
    if fmt == "fits":
        q["bands"] = "grz"
    return f"{SERVICE}/cutout.{fmt}?{urllib.parse.urlencode(q)}"


def viewer_url(ra, dec):
    return f"{SERVICE}?ra={ra:.6f}&dec={dec:.6f}&layer=ls-dr10&zoom=16"


def valid(body, fmt):
    """The service answers a bad request with HTTP 200 and a line of text."""
    return body.startswith(MAGIC[fmt])


def fetch(url, dest, fmt):
    """Download url to dest; None on success, else why not."""
    why = ""
    for wait in WAITS:
        time.sleep(wait)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read()
        except (urllib.error.URLError, OSError) as e:
            why = str(e)
            continue
        if not valid(body, fmt):
            why = f"not a {fmt}: {body[:80]!r}"
            continue
        with open(dest + ".part", "wb") as f:
            f.write(body)
        os.replace(dest + ".part", dest)
        return None
    return why


def have(path, fmt):
    if not os.path.exists(path):
        return False
    with open(path, "rb") as f:
        return valid(f.read(16), fmt)


def folder(cat, i):
    return cat.row(i)[cat.index["name"]].replace(" ", "_")


def jobs(cat, idx, a):
    """(entry, panel, url, path) for every cutout still to fetch."""
    ra, dec, r = cat.numeric("ra"), cat.numeric("dec"), cat.numeric("shape_r")
    out = []
    for i in idx:
        size, scale = field(r[i], a)
        panels = [(p, LAYERS[p], size, scale) for p in a.layers]
        if a.context:
            panels.append(("context", LAYERS["image"], size, 4 * scale))
        for panel, layer, s, sc in panels:
            path = os.path.join(a.out, folder(cat, i), f"{panel}.{a.format}")
            if not have(path, a.format):
                out.append((i, panel, cutout_url(layer, ra[i], dec[i], s, sc, a.format), path))
    return out


def gallery(cat, idx, cell, a):
    ra, dec, r = cat.numeric("ra"), cat.numeric("dec"), cat.numeric("shape_r")
    panels = list(a.layers) + (["context"] if a.context else [])
    rows = []
    for i in idx:
        row = cat.row(i)
        name = row[cat.index["name"]]
        size, scale = field(r[i], a)
        facts = " · ".join(f"{c} {html.escape(row[cat.index[c]])}" for c in SHOWN if c in cat.index)
        if i in cell:
            facts += f" · cell {html.escape(cell[i])}"
        imgs = []
        for p in panels:
            src = f"{folder(cat, i)}/{p}.{a.format}"
            fov = size * (4 * scale if p == "context" else scale)
            if a.format == "jpg":
                imgs.append(f'<figure><img loading="lazy" src="{src}" alt="{p}">'
                            f"<figcaption>{p} {fov:.0f}\"</figcaption></figure>")
            else:
                imgs.append(f'<figure><a href="{src}">{p}.fits</a><figcaption>{fov:.0f}"</figcaption></figure>')
        rows.append(f'<section><h2><a href="{viewer_url(ra[i], dec[i])}">{html.escape(name)}</a>'
                    f'<span>{facts}</span></h2><div class="row">{"".join(imgs)}</div></section>')
    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>DEBBIE images</title>
<style>
:root{{--bg:#fff;--fg:#1b1b1b;--dim:#666;--line:#ddd}}
@media (prefers-color-scheme: dark){{:root{{--bg:#121212;--fg:#e8e8e8;--dim:#999;--line:#333}}}}
body{{background:var(--bg);color:var(--fg);font:14px/1.4 system-ui,sans-serif;margin:0 auto;max-width:1100px;padding:16px}}
header p{{color:var(--dim)}}
section{{border-top:1px solid var(--line);padding:10px 0}}
h2{{font-size:15px;margin:0 0 6px;display:flex;flex-wrap:wrap;gap:4px 12px;align-items:baseline}}
h2 a{{color:inherit}} h2 span{{font-weight:400;color:var(--dim);font-size:13px}}
.row{{display:flex;flex-wrap:wrap;gap:8px}}
figure{{margin:0}} img{{width:220px;height:220px;display:block;background:#000}}
figcaption{{color:var(--dim);font-size:12px}}
</style></head><body>
<header><h1>DEBBIE images</h1>
<p>{len(idx):,} entries. Image, Tractor model and residual from Legacy Surveys DR10 at
{PIXSCALE}"/pixel. Each name links to the Legacy Surveys viewer.
Selection: <code>{html.escape(" ".join(sys.argv[1:]) or "(none)")}</code></p></header>
{"".join(rows)}
</body></html>
"""
    with open(os.path.join(a.out, "index.html"), "w", encoding="utf-8") as f:
        f.write(page)


def write_selection(cat, idx, cell, weight, a):
    with open(os.path.join(a.out, "selection.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["cell", "weight"] + cat.header)
        for i in idx:
            w.writerow([cell.get(i, ""), f"{weight.get(i, 1.0):.4f}"] + cat.row(i))


def select(cat, a):
    idx = list(range(len(cat)))
    numbers = [number_of(n) for n in a.name or []] + list(a.number or [])
    if numbers:
        idx = by_number(cat, idx, numbers)
    if a.cone:
        idx = cone(cat, idx, *a.cone)
    if a.part:
        idx = where(cat, idx, f"part=={a.part}")
    for expr in a.where or []:
        idx = where(cat, idx, expr)
    print(f"selected {len(idx):,} of {len(cat):,}")
    cell, weight = {}, {}
    rng = random.Random(a.seed)
    if a.stratify:
        if not a.sample:
            raise SystemExit("--stratify needs --sample N")
        idx, cell, weight = stratified(cat, idx, a.stratify, a.sample, rng)
    elif a.sample and a.sample < len(idx):
        n = len(idx)
        idx = rng.sample(idx, a.sample)
        weight = {i: n / a.sample for i in idx}
    if a.sort:
        idx = order(cat, idx, a.sort, a.desc)
    if a.limit is not None:
        idx = idx[:a.limit]
    return idx, cell, weight


def parser():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog="\n\n".join(__doc__.split("\n\n")[1:]))
    ap.add_argument("--catalogue", default=CATALOGUE, help="the DEBBIE CSV, gzipped or not")
    ap.add_argument("--name", nargs="+", metavar="NAME")
    ap.add_argument("--number", nargs="+", type=int, metavar="N")
    ap.add_argument("--cone", nargs=3, type=float, metavar=("RA", "DEC", "ARCMIN"))
    ap.add_argument("--part", choices=["original", "recovered", "extension"])
    ap.add_argument("--where", nargs="+", metavar="EXPR")
    ap.add_argument("--sample", type=int, metavar="N")
    ap.add_argument("--stratify", nargs="+", metavar="COL:K")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--sort", metavar="COL")
    ap.add_argument("--desc", action="store_true", help="with --sort: largest first")
    ap.add_argument("--limit", type=int, metavar="N")
    ap.add_argument("--layers", nargs="+", choices=list(LAYERS), default=list(LAYERS))
    ap.add_argument("--context", action="store_true", help="add an image four times wider")
    ap.add_argument("--format", choices=["jpg", "fits"], default="jpg")
    ap.add_argument("--fov-min", type=float, default=40.0, metavar="ARCSEC")
    ap.add_argument("--fov-reff", type=float, default=8.0, metavar="X")
    ap.add_argument("--out", default="debbie_images")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--max-objects", type=int, default=1000)
    ap.add_argument("--dry-run", action="store_true", help="write selection.csv only")
    ap.add_argument("--columns", action="store_true", help="list the catalogue's columns")
    return ap


def main(argv=None):
    a = parser().parse_args(argv)
    if not 1 <= a.workers <= 8:
        raise SystemExit("--workers must be 1 to 8; the cutout service is shared")
    cat = Catalogue(a.catalogue)
    if a.columns:
        print("\n".join(cat.header))
        return
    idx, cell, weight = select(cat, a)
    if not idx:
        raise SystemExit("nothing selected")
    if len(idx) > a.max_objects and not a.dry_run:
        raise SystemExit(f"{len(idx):,} entries is more than --max-objects {a.max_objects:,}; "
                         "use --sample or --limit, or raise --max-objects")
    os.makedirs(a.out, exist_ok=True)
    write_selection(cat, idx, cell, weight, a)
    if a.dry_run:
        print(f"wrote {a.out}/selection.csv ({len(idx):,} entries); nothing downloaded")
        return
    todo = jobs(cat, idx, a)
    print(f"{len(idx):,} entries, {len(todo):,} cutouts to fetch into {a.out}/")
    failed = []
    with ThreadPoolExecutor(a.workers) as pool:
        futures = {}
        for i, panel, url, path in todo:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            futures[pool.submit(fetch, url, path, a.format)] = (i, panel, url)
        for k, fut in enumerate(as_completed(futures), 1):
            why = fut.result()
            if why:
                i, panel, url = futures[fut]
                failed.append((cat.row(i)[cat.index["name"]], panel, url, why))
            if k % 50 == 0 or k == len(todo):
                print(f"  {k:,}/{len(todo):,} fetched, {len(failed)} failed", flush=True)
    if failed:
        with open(os.path.join(a.out, "failures.csv"), "w", newline="", encoding="utf-8") as f:
            csv.writer(f).writerows([("name", "panel", "url", "why")] + failed)
        print(f"{len(failed)} cutouts failed: {a.out}/failures.csv; run again to retry them")
    gallery(cat, idx, cell, a)
    print(f"gallery: {os.path.abspath(os.path.join(a.out, 'index.html'))}")


if __name__ == "__main__":
    main()
