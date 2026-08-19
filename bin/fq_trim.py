#!/usr/bin/env python3
"""Quality/adapter trimming for FASTQ. Standard library only.

Trailing-quality trim (running-sum, cutadapt-style) then adapter match then
length filter. Deterministic: no randomness, so output is byte-identical
across runs on the same input.
"""
import argparse
import gzip
import sys

from fq_qc import iter_fastq, open_maybe_gzip, phred_scores


def trim_trailing_quality(seq, qual, cutoff):
    """cutadapt-style running-sum trim from the 3' end.

    Walks right-to-left accumulating (cutoff - q); the cut point is where that
    running sum is maximal. Tolerates a single low base inside a good tail
    rather than truncating at the first dip.
    """
    qs = phred_scores(qual)
    running = 0
    max_running = 0
    cut = len(seq)
    for i in range(len(qs) - 1, -1, -1):
        running += cutoff - qs[i]
        if running < 0:
            break
        if running > max_running:
            max_running = running
            cut = i
    return seq[:cut], qual[:cut]


def find_adapter(seq, adapter, min_overlap=3, max_mismatch_rate=0.15):
    """Return the 5'-most index where `adapter` matches, else -1.

    Also matches partial adapter at the read end (3' overhang), which is where
    read-through adapter actually appears.
    """
    n, m = len(seq), len(adapter)
    for start in range(n):
        overlap = min(m, n - start)
        if overlap < min_overlap:
            break
        mismatches = sum(
            1 for k in range(overlap) if seq[start + k] != adapter[k]
        )
        if mismatches <= int(overlap * max_mismatch_rate):
            return start
    return -1


def main(argv=None):
    ap = argparse.ArgumentParser(description="Dependency-free FASTQ trimmer")
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out-fastq", required=True)
    ap.add_argument("--out-log", required=True)
    ap.add_argument("--quality-cutoff", type=int, default=20)
    ap.add_argument("--min-length", type=int, default=25)
    ap.add_argument("--adapter", default="AGATCGGAAGAGC")
    args = ap.parse_args(argv)

    stats = {
        "sample": args.sample,
        "reads_in": 0,
        "reads_out": 0,
        "reads_dropped_too_short": 0,
        "reads_adapter_trimmed": 0,
        "reads_quality_trimmed": 0,
        "bases_in": 0,
        "bases_out": 0,
    }

    opener = gzip.open if args.out_fastq.endswith(".gz") else open
    with open_maybe_gzip(args.fastq) as fin, opener(args.out_fastq, "wt") as fout:
        for header, seq, qual in iter_fastq(fin):
            stats["reads_in"] += 1
            stats["bases_in"] += len(seq)

            idx = find_adapter(seq, args.adapter)
            if idx >= 0:
                seq, qual = seq[:idx], qual[:idx]
                stats["reads_adapter_trimmed"] += 1

            pre = len(seq)
            seq, qual = trim_trailing_quality(seq, qual, args.quality_cutoff)
            if len(seq) < pre:
                stats["reads_quality_trimmed"] += 1

            if len(seq) < args.min_length:
                stats["reads_dropped_too_short"] += 1
                continue

            stats["reads_out"] += 1
            stats["bases_out"] += len(seq)
            fout.write(f"{header}\n{seq}\n+\n{qual}\n")

    cols = list(stats.keys())
    with open(args.out_log, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        fh.write("\t".join(str(stats[c]) for c in cols) + "\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
