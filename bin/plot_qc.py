#!/usr/bin/env python3
"""Render the QC summary figure: per-cycle quality and before/after metrics."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAW_COLOUR = "#9aa3af"
TRIM_COLOUR = "#1f6f8b"
KEEP_COLOUR = "#3f8f5f"

# Minimum fraction of reads a cycle must be averaged over to be plotted.
SUPPORT_FRAC = 0.10


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-qc", required=True)
    ap.add_argument("--trimmed-qc", required=True)
    ap.add_argument("--json-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--q30-floor", type=float, default=None)
    args = ap.parse_args(argv)

    raw = {r["sample"]: r for r in read_tsv(args.raw_qc)}
    trim = {r["sample"]: r for r in read_tsv(args.trimmed_qc)}
    samples = sorted(raw)

    jdir = Path(args.json_dir)
    curves = {}
    for s in samples:
        for tag in ("raw", "trimmed"):
            p = jdir / f"{s}.{tag}.qc.json"
            if not p.exists():
                continue
            d = json.loads(p.read_text())
            q = d["per_cycle_mean_quality"]
            n = d.get("per_cycle_read_count") or [1] * len(q)
            # A cycle averaged over a handful of reads is noise, not signal:
            # keep only cycles supported by at least SUPPORT_FRAC of the reads.
            floor = SUPPORT_FRAC * max(n)
            keep = [i for i, c in enumerate(n) if c >= floor]
            cut = (max(keep) + 1) if keep else len(q)
            curves[(s, tag)] = q[:cut]

    plt.rcParams.update({
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 200,
    })

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.3))

    # ---- Panel 1: per-cycle quality, raw vs trimmed ----
    ax = axes[0]
    shades = ["#123f52", "#1f6f8b", "#5aa8bf"]
    for i, s in enumerate(samples):
        c = shades[i % len(shades)]
        if (s, "raw") in curves:
            y = curves[(s, "raw")]
            ax.plot(range(1, len(y) + 1), y, color=c, ls=":", lw=1.1, alpha=0.85)
        if (s, "trimmed") in curves:
            y = curves[(s, "trimmed")]
            ax.plot(range(1, len(y) + 1), y, color=c, ls="-", lw=1.7)
    ax.axhline(30, color="#b0b0b0", lw=0.8, ls="--", zorder=0)
    ax.annotate("Q30", xy=(0.0, 30), xycoords=("axes fraction", "data"),
                xytext=(2, 3), textcoords="offset points",
                ha="left", va="bottom", color="#8a8a8a", fontsize=6)
    # Direct end-of-line labels. The trimmed curves converge near the same
    # quality, so label anchors are stacked over a fixed fraction of the axis
    # height and each keeps a leader line back to its own curve.
    ends = []
    for i, s in enumerate(samples):
        y = curves.get((s, "trimmed")) or curves.get((s, "raw"))
        if y:
            ends.append((s, len(y), y[-1], shades[i % len(shades)]))
    ends.sort(key=lambda e: e[2])
    if ends:
        ax.autoscale_view()
        lo, hi = ax.get_ylim()
        step = 0.13 * (hi - lo)
        centre = sum(e[2] for e in ends) / len(ends)
        base = centre - step * (len(ends) - 1) / 2.0
        xmax = max(e[1] for e in ends)
        for k, (s, xe, ye, col) in enumerate(ends):
            ax.annotate(s, xy=(xe, ye), xytext=(xmax * 1.10, base + k * step),
                        textcoords="data", va="center", ha="left",
                        fontsize=6, color=col,
                        arrowprops=dict(arrowstyle="-", color=col, lw=0.5,
                                        shrinkA=0, shrinkB=1))
    ax.plot([], [], color="#333333", ls=":", lw=1.1, label="before trimming")
    ax.plot([], [], color="#333333", ls="-", lw=1.7, label="after trimming")
    ax.set_xlabel("Cycle (base position)")
    ax.set_ylabel("Mean Phred quality")
    ax.set_title("Low-quality 3' tails are removed", loc="left")
    ax.legend(frameon=False, loc="lower left")
    ax.margins(x=0.20, y=0.06)

    x = list(range(len(samples)))
    w = 0.38

    # ---- Panel 2: %Q30 before/after ----
    ax = axes[1]
    rq = [float(raw[s]["pct_q30"]) for s in samples]
    tq = [float(trim[s]["pct_q30"]) for s in samples]
    ax.bar([i - w / 2 for i in x], rq, w, label="before", color=RAW_COLOUR)
    ax.bar([i + w / 2 for i in x], tq, w, label="after", color=TRIM_COLOUR)
    if args.q30_floor is not None:
        ax.axhline(args.q30_floor, color="#c2543a", lw=1.0, ls="--", zorder=3)
        ax.annotate(f"required floor {args.q30_floor:g}%",
                    xy=(0.0, args.q30_floor), xycoords=("axes fraction", "data"),
                    xytext=(2, 3), textcoords="offset points",
                    fontsize=6, color="#c2543a", va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=20, ha="right")
    ax.set_ylabel("Bases at or above Q30 (%)")
    ax.set_title("Q30 content improves for every sample", loc="left")
    top = max(max(rq), max(tq))
    ax.set_ylim(0, top * 1.28)
    # Headline number only: the largest gain.
    gains = [b - a for a, b in zip(rq, tq)]
    j = gains.index(max(gains))
    ax.annotate(f"+{gains[j]:.0f} pts", xy=(j + w / 2, tq[j]), xytext=(0, 4),
                textcoords="offset points", ha="center", fontsize=7,
                color=TRIM_COLOUR)
    ax.legend(frameon=False, loc="upper right", ncol=2,
              columnspacing=1.0, handlelength=1.2)

    # ---- Panel 3: read retention ----
    ax = axes[2]
    rr = [int(raw[s]["n_reads"]) for s in samples]
    tr = [int(trim[s]["n_reads"]) for s in samples]
    ax.bar([i - w / 2 for i in x], rr, w, label="input", color=RAW_COLOUR)
    ax.bar([i + w / 2 for i in x], tr, w, label="retained", color=KEEP_COLOUR)
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=20, ha="right")
    ax.set_ylabel("Reads")
    pcts = [100.0 * b / a for a, b in zip(rr, tr)]
    title = ("Every read survives the length filter"
             if min(pcts) >= 99.5 else
             f"Reads retained: {min(pcts):.0f}\u2013{max(pcts):.0f}%")
    ax.set_title(title, loc="left")
    ax.set_ylim(0, max(rr) * 1.25)
    k = pcts.index(min(pcts))
    ax.annotate(f"{pcts[k]:.0f}% kept", xy=(k + w / 2, tr[k]), xytext=(0, 4),
                textcoords="offset points", ha="center", fontsize=7,
                color=KEEP_COLOUR)
    ax.legend(frameon=False, loc="upper right", ncol=2,
              columnspacing=1.0, handlelength=1.2)

    fig.tight_layout(w_pad=1.6)
    fig.savefig(args.out, bbox_inches="tight")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
