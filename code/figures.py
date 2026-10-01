"""Build the paper figures.

Fig 1  the radial representation: an object, its annuli, its profile
Fig 2  benchmark AUC, ours vs the published CNN, with the ablation
Fig 3  sky distribution and the selection function
Fig 4  montage: selected objects above, a random draw below

Everything is drawn from files in the repository and from the v1.2 catalogue,
so the figures regenerate from the same numbers the text quotes.

Design rules applied throughout, and why:

  * ONE accent hue. Fig 2's story is a single number -- the gain from the
    representation -- so the winning bar carries the accent and everything
    else is gray. Three arbitrary hues would spend the colour channel on
    identity nobody needs.
  * The selection-function heatmap is a ONE-HUE ramp, light to dark. It was
    `magma`, which is multi-hue and therefore not monotonic in perceived
    lightness once printed in greyscale; a single-hue ramp survives a mono
    printer, which is still how referees read.
  * Hairline axes, no top/right spines, no dashed gridlines. Dashes read as
    "threshold" when they are only a grid.
  * Labels are placed selectively -- the extremes and the winner -- never one
    on every mark.

Fig 3a keeps its dark treatment: at this density a light ground turns the
footprint into undifferentiated ink, and the galactic-plane exclusion stops
being legible, which is the whole point of the panel.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, "/home/ben/star")
sys.path.insert(0, "/home/ben/star/src")
sys.path.insert(0, "/home/ben/legacy")

LEGACY = "/home/ben/legacy"
FIGS = "/home/ben/star/paper/figures"
CAT = "/home/ben/star/catalogue_debbie/v1.2_main.npy"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"            # primary text
INK2 = "#52514e"           # secondary text
RULE = "#dcdbd6"           # hairline grid / axis
ACCENT = "#2a78d6"         # the one hue that carries meaning
MUTED = "#8a8985"          # everything not being emphasised

# one-hue sequential ramp, light -> dark (greyscale-safe: lightness monotonic)
SEQ = LinearSegmentedColormap.from_list("seqblue", [
    "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b",
])

plt.rcParams.update({
    "font.size": 8.5,
    "axes.titlesize": 9,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "figure.dpi": 400,
    "savefig.dpi": 400,
    "savefig.bbox": "tight",
    "savefig.facecolor": SURFACE,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "axes.grid": False,
    "axes.edgecolor": RULE,
    "axes.linewidth": 0.6,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "font.family": "DejaVu Sans",
})


def _clean(ax, left=True, bottom=True):
    """Hairline frame: drop the top and right rules, keep the data ones."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)
    ax.tick_params(length=3, width=0.6)


def img_path(ra, dec, arcsec=75):
    return (f"{LEGACY}/candidate_images/ra{int(ra):03d}/"
            f"{ra:.5f}{dec:+.5f}_{arcsec}.jpg")


def fig1_representation():
    """An object, the annuli laid over it, and the profile that results."""
    from PIL import Image
    # The example must be an object the catalogue still contains. The
    # previous pick (299.4597-30.6674) was in the low-latitude extension AND
    # was later removed by the DESI screen, so the figure illustrating the
    # method showed an object no longer in the catalogue.
    sl = np.load("/home/ben/star/catalogue_debbie/DEBBIE_SHORTLIST.npy")
    # Chosen by measurement, not by eye: lowest flux-centroid offset from the
    # frame centre (0.2%) with the highest core contrast (10.2 sigma) among
    # the 40 largest shortlist objects. An off-centre example makes panel (b)
    # look like the annuli are misaligned with the object, which is exactly
    # the wrong thing for the figure defining the method to imply.
    want = (61.80906, 5.02184)
    hit = np.flatnonzero((np.abs(sl["ra"] - want[0]) < 1e-3)
                         & (np.abs(sl["dec"] - want[1]) < 1e-3))
    o = sl[hit[0]] if len(hit) else sl[np.argsort(-sl["shape_r"])][0]
    p = img_path(float(o["ra"]), float(o["dec"]))
    if not os.path.exists(p):
        print("fig1: no image"); return
    im = np.asarray(Image.open(p).convert("RGB"))
    g = np.asarray(Image.open(p).convert("L"), dtype=float)
    n = g.shape[0]
    c = n / 2.0

    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
    ax[0].imshow(im)
    ax[0].set_title("(a) cutout", loc="left", color=INK)
    ax[0].axis("off")

    ax[1].imshow(im)
    r_edges = np.geomspace(2, n / 2, 12)
    for r in r_edges:
        ax[1].add_patch(plt.Circle((c, c), r, fill=False, color=ACCENT,
                                   lw=0.45, alpha=0.85))
    for k in range(16):
        a = 2 * np.pi * k / 16
        ax[1].plot([c, c + n / 2 * np.cos(a)], [c, c + n / 2 * np.sin(a)],
                   color=ACCENT, lw=0.28, alpha=0.55)
    ax[1].set_xlim(0, n); ax[1].set_ylim(n, 0)
    ax[1].set_title(r"(b) annuli $\times$ position angle", loc="left", color=INK)
    ax[1].axis("off")

    yy, xx = np.mgrid[0:n, 0:n]
    rr = np.hypot(yy - c, xx - c)
    prof, err, rc = [], [], []
    for i in range(len(r_edges) - 1):
        m = (rr >= r_edges[i]) & (rr < r_edges[i + 1])
        if m.sum() < 5:
            continue
        v = g[m]
        prof.append(np.median(v)); err.append(np.std(v) / np.sqrt(m.sum()))
        rc.append(0.5 * (r_edges[i] + r_edges[i + 1]))
    bg = np.median(g[rr > n * 0.45])
    prof = np.array(prof) - bg
    ax[2].axhline(0, color=RULE, lw=0.6, zorder=1)
    ax[2].errorbar(rc, prof, yerr=err, marker="o", ms=3.2, lw=1.4,
                   color=ACCENT, ecolor=ACCENT, elinewidth=0.8, capsize=0,
                   markeredgecolor=SURFACE, markeredgewidth=0.6, zorder=3)
    ax[2].set_xscale("log")
    ax[2].set_xlabel("radius (pixels)")
    ax[2].set_ylabel("surface brightness $-$ background")
    ax[2].set_title("(c) radial profile", loc="left", color=INK)
    _clean(ax[2])
    fig.tight_layout(w_pad=1.6)
    fig.savefig(f"{FIGS}/fig1_representation.png")
    plt.close(fig)
    print("fig1 done")


def fig2_benchmark():
    """The benchmark result and the ablation that explains it.

    A dot-and-interval plot, not bars. Every AUC here sits between 0.90 and
    0.98, so a bar chart would have to start its axis at 0.88 -- and a bar
    whose baseline is not zero overstates the ratio between values by however
    much was cut off. A dot carries no baseline claim, so the axis can span
    the range that matters honestly.

    Emphasis, not a palette: the representation is the finding, so it alone
    carries the accent hue. The published CNN is drawn hollow because it was
    not rerun here.
    """
    b = json.load(open("/home/ben/star/output/benchmark_deepshadows.json"))
    rows = [
        ("six catalogue features", "", b["six_features"]["auc"],
         b["six_features"]["ci"], MUTED, False),
        ("DeepShadows CNN", "as published", 0.9740, (0.9700, 0.9740),
         MUTED, True),
        ("radial representation", "", b["radial"]["auc"],
         b["radial"]["ci"], ACCENT, False),
    ]
    fig, ax = plt.subplots(figsize=(5.6, 2.15))
    y = np.arange(len(rows))
    for i, (_, _, v, ci, col, published) in enumerate(rows):
        ax.plot([ci[0], ci[1]], [i, i], color=col, lw=1.6, zorder=2,
                solid_capstyle="round", alpha=0.55)
        ax.plot([v], [i], marker="o", ms=8, zorder=3,
                color=SURFACE if published else col,
                markeredgecolor=col, markeredgewidth=1.6)
        ax.text(1.001, i, f"{v:.4f}", va="center", ha="left", fontsize=8,
                color=INK, family="DejaVu Sans Mono",
                fontweight="bold" if col == ACCENT else "normal",
                transform=ax.get_yaxis_transform(), clip_on=False, zorder=4)
    labels = [name if not sub else f"{name}  ({sub})"
              for name, sub, *_ in rows]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, color=INK, fontsize=8)
    for t, (_, _, _, _, col, _) in zip(ax.get_yticklabels(), rows):
        if col == ACCENT:
            t.set_fontweight("bold")
    ax.set_ylim(-0.45, len(rows) - 0.15)
    ax.set_xlim(0.895, 0.995)
    ax.set_xticks([0.90, 0.92, 0.94, 0.96, 0.98])
    ax.set_xlabel("AUC on the DeepShadows benchmark (hard negatives)")
    ax.set_title("Identical classifier and pipeline except where noted",
                 loc="left", fontsize=7.5, color=INK2, pad=2)
    ax.xaxis.grid(True, color=RULE, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    _clean(ax, left=False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    # The span is the whole point of the figure, so connect it to the two
    # dots it actually spans -- an arrow floating above unlabelled rows
    # invites the reader to pair it with the wrong ones.
    lo, hi = b["six_features"]["auc"], b["radial"]["auc"]
    yb = 2.30
    for x, y0 in ((lo, 0), (hi, 2)):
        ax.plot([x, x], [y0, yb], color=RULE, lw=0.6, zorder=1)
    ax.annotate("", xy=(hi, yb), xytext=(lo, yb),
                arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.8,
                                shrinkA=0, shrinkB=0))
    ax.text(0.5 * (lo + hi), yb + 0.07,
            f"+{hi - lo:.3f} from the representation alone",
            ha="center", va="bottom", fontsize=7.5, color=INK2)
    fig.tight_layout()
    fig.subplots_adjust(right=0.86)
    fig.savefig(f"{FIGS}/fig2_benchmark.png")
    plt.close(fig)
    print("fig2 done")


def fig3_sky_and_selection():
    """Where the catalogue is, and what fraction of galaxies the cuts keep."""
    cat = np.load(CAT)
    comp = json.load(open(f"{LEGACY}/completeness.json"))

    SKY_INK = "#0d1017"
    SKY = "#7fc8f8"
    SKY_GRID = "#5a6472"

    fig = plt.figure(figsize=(7.2, 3.1))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.22,
                          top=0.86, bottom=0.16, left=0.04, right=0.97)

    ax1 = fig.add_subplot(gs[0, 0], projection="mollweide")
    ax1.set_facecolor(SKY_INK)
    ra = np.radians(((cat["ra"] + 180) % 360) - 180)
    dec = np.radians(cat["dec"])
    # Sub-pixel markers at low alpha: density then builds where the sky is
    # genuinely crowded, and single objects stay visible at the sparse edges.
    ax1.scatter(ra, dec, s=0.05, c=SKY, alpha=0.30, linewidths=0,
                rasterized=True)
    ax1.grid(True, lw=0.35, alpha=0.5, color=SKY_GRID)
    ax1.set_xticklabels([])
    ax1.tick_params(colors=SKY_GRID, labelsize=6.5)
    for sp in ax1.spines.values():
        sp.set_color(SKY_GRID); sp.set_linewidth(0.7)
    # Titles are placed in figure coordinates so the two panels align; a
    # mollweide axes title floats with the projection and will not.
    fig.text(0.055, 0.955, f"(a) {len(cat):,} candidates, DR10 south, "
             r"$|b| > 30^\circ$", fontsize=8.5, color=INK, ha="left",
             va="top")

    ax2 = fig.add_subplot(gs[0, 1])
    grid = comp["grid"]
    rs = sorted({k.split("_")[0] for k in grid},
                key=lambda s: float(s.split("-")[0]))
    # Only the mu columns the selection actually targets. Carrying tick
    # labels for 22-23 or 26.5-30 advertises columns the cuts exclude by
    # construction, which reads as missing data rather than as out of scope.
    ms = sorted({k.split("_")[1] for k in grid
                 if 24.0 <= float(k.split("_")[1].split("-")[0]) < 26.5},
                key=lambda s: float(s.split("-")[0]))
    M = np.full((len(rs), len(ms)), np.nan)
    for i, r in enumerate(rs):
        for j, m in enumerate(ms):
            k = f"{r}_{m}"
            if k in grid:
                M[i, j] = grid[k]["frac"] * 100
    Mm = M
    im = ax2.imshow(Mm, aspect="auto", origin="lower", cmap=SEQ,
                    vmin=0, vmax=100)
    ax2.set_xticks(range(len(ms)))
    ax2.set_xticklabels(ms, fontsize=6.5)
    ax2.set_yticks(range(len(rs)))
    ax2.set_yticklabels(rs, fontsize=6.5)
    ax2.set_xlabel(r"$\mu_{\rm eff}$ (mag arcsec$^{-2}$)")
    ax2.set_ylabel(r"$r_{\rm eff}$ (arcsec)")
    fig.text(0.575, 0.955, "(b) selection fraction, photometric cuts only",
             fontsize=8.5, color=INK, ha="left", va="top")
    fig.text(0.575, 0.905, "blank cells hold fewer than 30 reference galaxies",
             fontsize=6.8, color=INK2, ha="left", va="top")
    for i in range(len(rs)):
        for j in range(len(ms)):
            if np.isfinite(Mm[i, j]):
                ax2.text(j, i, f"{Mm[i, j]:.0f}", ha="center", va="center",
                         fontsize=6.2,
                         color="#ffffff" if Mm[i, j] > 55 else INK)
    for sp in ax2.spines.values():
        sp.set_visible(False)
    ax2.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax2, fraction=0.045, pad=0.03)
    cb.set_label("% retained", fontsize=7)
    cb.ax.tick_params(labelsize=6.5, length=2, width=0.5)
    cb.outline.set_visible(False)
    fig.savefig(f"{FIGS}/fig3_sky_selection.png")
    plt.close(fig)
    print("fig3 done")


def fig4_montage():
    """Selected objects above the rule, a random draw below it.

    A montage of only the best objects tells the reader what the pipeline can
    find, never what it typically returns. Showing both, labelled, is the
    honest version and costs one extra row.
    """
    from PIL import Image
    sl = np.load("/home/ben/star/catalogue_debbie/DEBBIE_SHORTLIST.npy")
    sl = sl[np.argsort(-sl["shape_r"])]
    cat = np.load(CAT)

    sel, i = [], 0
    while len(sel) < 8 and i < len(sl):
        o = sl[i]; i += 1
        if os.path.exists(img_path(float(o["ra"]), float(o["dec"]))):
            sel.append(o)

    rng = np.random.default_rng(20260905)
    rnd, tries = [], 0
    while len(rnd) < 8 and tries < 4000:
        o = cat[rng.integers(len(cat))]; tries += 1
        if os.path.exists(img_path(float(o["ra"]), float(o["dec"]))):
            rnd.append(o)

    fig = plt.figure(figsize=(7.2, 7.9))
    outer = fig.add_gridspec(2, 1, hspace=0.13, top=0.935, bottom=0.015,
                             left=0.03, right=0.985)
    blocks = [("Selected examples: the largest objects that clear every screen",
               sel, 0),
              ("Random draw from the catalogue, at the same scale", rnd, 1)]
    for header, objs, bi in blocks:
        gs = outer[bi].subgridspec(2, 4, hspace=0.18, wspace=0.05)
        pos = outer[bi].get_position(fig)
        fig.text(0.03, pos.y1 + 0.012, header, fontsize=8, color=INK,
                 ha="left", va="bottom")
        for k, o in enumerate(objs):
            ax = fig.add_subplot(gs[k // 4, k % 4])
            ax.imshow(np.asarray(Image.open(
                img_path(float(o["ra"]), float(o["dec"]))).convert("RGB")))
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_color(RULE); sp.set_linewidth(0.5)
            ax.set_title(f"{o['ra']:.4f}{o['dec']:+.4f}", fontsize=6.0,
                         color=INK2, pad=2.5)
            # size label on a translucent plate: several cutouts have a bright
            # star in the corner that swallows plain white text
            ax.text(0.04, 0.04, f"{o['shape_r']:.1f}\u2033",
                    transform=ax.transAxes, fontsize=6.2, color="#ffffff",
                    ha="left", va="bottom",
                    bbox=dict(boxstyle="round,pad=0.18", facecolor="#000000",
                              edgecolor="none", alpha=0.55))
    fig.savefig(f"{FIGS}/fig4_montage.png")
    plt.close(fig)
    print(f"fig4 done ({len(sel)} selected + {len(rnd)} random)")


if __name__ == "__main__":
    os.makedirs(FIGS, exist_ok=True)
    fig1_representation()
    fig2_benchmark()
    fig3_sky_and_selection()
    fig4_montage()
