#!/usr/bin/env python3
"""Assemble the run's QC tables and validation status into a markdown report."""
import argparse
import csv
from pathlib import Path


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def md_table(rows, cols):
    out = ["| " + " | ".join(cols) + " |",
           "|" + "|".join(["---"] * len(cols)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-qc", required=True)
    ap.add_argument("--trimmed-qc", required=True)
    ap.add_argument("--trim-log", required=True)
    ap.add_argument("--validation", required=True)
    ap.add_argument("--figure", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    raw = read_tsv(args.raw_qc)
    trim = read_tsv(args.trimmed_qc)
    log = read_tsv(args.trim_log)
    validation = Path(args.validation).read_text().strip()

    cols = ["sample", "n_reads", "mean_read_length", "mean_quality",
            "pct_q20", "pct_q30", "gc_percent"]
    logcols = ["sample", "reads_in", "reads_out", "reads_adapter_trimmed",
               "reads_quality_trimmed", "reads_dropped_too_short"]

    body = f"""# FASTQ QC Report

## Validation status

```
{validation}
```

The validation step is a hard gate. The workflow exits non-zero if trimming
increases read or base counts, if mean quality or %Q30 decreases, if trimmed
%Q30 falls below the configured floor, or if the trim logs fail to reconcile
with the QC read counts. A report is only produced when every invariant holds.

## Summary figure

![QC summary]({Path(args.figure).name})

## Raw FASTQ metrics

{md_table(raw, cols)}

## Trimmed FASTQ metrics

{md_table(trim, cols)}

## Trimming statistics

{md_table(log, logcols)}

## Methods

Reads were assessed for per-base quality, GC content, read-length distribution
and ambiguous-base content using a dependency-free Python implementation
(`bin/fq_qc.py`). Adapter sequence was removed by mismatch-tolerant prefix
matching (default 15 percent mismatch rate, minimum 3 bp overlap), followed by
3'-end quality trimming using a cutadapt-style running-sum algorithm at the
configured Phred cutoff, and a minimum-length filter (`bin/fq_trim.py`).
Trimmed reads were re-assessed with the identical QC implementation so that
before/after metrics are directly comparable. All steps are deterministic:
repeated runs on identical input produce byte-identical output.
"""
    Path(args.out).write_text(body)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
