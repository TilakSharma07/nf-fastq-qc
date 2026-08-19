#!/usr/bin/env python3
"""Per-sample FASTQ QC metrics. Standard library only — no third-party deps.

Emits a single-row TSV so downstream processes can concatenate sample metrics
without a parsing dependency.
"""
import argparse
import gzip
import json
import sys
from collections import Counter

PHRED_OFFSET = 33


def open_maybe_gzip(path):
    """Open plain or gzipped FASTQ transparently, in text mode."""
    if str(path).endswith((".gz", ".gzip")):
        return gzip.open(path, "rt")
    return open(path, "rt")


def iter_fastq(handle):
    """Yield (header, sequence, quality) tuples, validating record structure.

    Raises ValueError on malformed records rather than silently skipping them,
    so a truncated file fails the pipeline instead of reaching the report.
    """
    while True:
        header = handle.readline()
        if not header:
            return
        header = header.rstrip("\n")
        seq = handle.readline().rstrip("\n")
        plus = handle.readline().rstrip("\n")
        qual = handle.readline().rstrip("\n")

        if not header.startswith("@"):
            raise ValueError(f"record header does not start with '@': {header[:40]!r}")
        if not plus.startswith("+"):
            raise ValueError(f"missing '+' separator after sequence for {header[:40]!r}")
        if len(seq) != len(qual):
            raise ValueError(
                f"sequence/quality length mismatch for {header[:40]!r}: "
                f"{len(seq)} vs {len(qual)}"
            )
        yield header, seq, qual


def phred_scores(qual):
    return [ord(c) - PHRED_OFFSET for c in qual]


def compute_metrics(path, sample):
    n_reads = 0
    total_bases = 0
    total_q = 0
    q20_bases = 0
    q30_bases = 0
    base_counts = Counter()
    lengths = []
    per_cycle_q_sum = []
    per_cycle_n = []

    with open_maybe_gzip(path) as fh:
        for _header, seq, qual in iter_fastq(fh):
            n_reads += 1
            L = len(seq)
            lengths.append(L)
            total_bases += L
            base_counts.update(seq.upper())

            qs = phred_scores(qual)
            total_q += sum(qs)
            q20_bases += sum(1 for q in qs if q >= 20)
            q30_bases += sum(1 for q in qs if q >= 30)

            if L > len(per_cycle_q_sum):
                pad = L - len(per_cycle_q_sum)
                per_cycle_q_sum.extend([0] * pad)
                per_cycle_n.extend([0] * pad)
            for i, q in enumerate(qs):
                per_cycle_q_sum[i] += q
                per_cycle_n[i] += 1

    if n_reads == 0:
        raise ValueError(f"{path}: no FASTQ records found")

    gc = base_counts["G"] + base_counts["C"]
    at = base_counts["A"] + base_counts["T"]
    gc_pct = 100.0 * gc / (gc + at) if (gc + at) else 0.0

    per_cycle_mean_q = [
        round(s / n, 4) for s, n in zip(per_cycle_q_sum, per_cycle_n) if n
    ]

    return {
        "sample": sample,
        "n_reads": n_reads,
        "total_bases": total_bases,
        "mean_read_length": round(total_bases / n_reads, 3),
        "min_read_length": min(lengths),
        "max_read_length": max(lengths),
        "mean_quality": round(total_q / total_bases, 4),
        "pct_q20": round(100.0 * q20_bases / total_bases, 4),
        "pct_q30": round(100.0 * q30_bases / total_bases, 4),
        "gc_percent": round(gc_pct, 4),
        "n_bases_ambiguous": base_counts["N"],
        "per_cycle_mean_quality": per_cycle_mean_q,
        # How many reads contributed to each cycle. After trimming, read lengths
        # vary, so late cycles are averages over few reads -- consumers should
        # suppress cycles with thin support rather than plotting noise.
        "per_cycle_read_count": [n for n in per_cycle_n if n],
    }


TSV_COLUMNS = [
    "sample", "n_reads", "total_bases", "mean_read_length",
    "min_read_length", "max_read_length", "mean_quality",
    "pct_q20", "pct_q30", "gc_percent", "n_bases_ambiguous",
]


def main(argv=None):
    ap = argparse.ArgumentParser(description="Dependency-free FASTQ QC metrics")
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-json", required=True)
    args = ap.parse_args(argv)

    m = compute_metrics(args.fastq, args.sample)

    with open(args.out_tsv, "w") as fh:
        fh.write("\t".join(TSV_COLUMNS) + "\n")
        fh.write("\t".join(str(m[c]) for c in TSV_COLUMNS) + "\n")

    with open(args.out_json, "w") as fh:
        json.dump(m, fh, indent=2)

    return 0


if __name__ == "__main__":
    sys.exit(main())
