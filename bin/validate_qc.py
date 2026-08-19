#!/usr/bin/env python3
"""Assert QC invariants; exit non-zero on violation.

This is the pipeline's positive control. Trimming must improve quality and must
not invent reads or bases. If any invariant fails, the run stops here rather
than producing a report that looks fine.
"""
import argparse
import csv
import sys


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def check(rows_raw, rows_trim, rows_log, min_q30_after=70.0):
    failures = []
    raw = {r["sample"]: r for r in rows_raw}
    trim = {r["sample"]: r for r in rows_trim}
    log = {r["sample"]: r for r in rows_log}

    if set(raw) != set(trim):
        failures.append(
            f"sample sets differ between raw and trimmed QC: "
            f"{sorted(set(raw) ^ set(trim))}"
        )

    for s in sorted(set(raw) & set(trim)):
        r, t, lg = raw[s], trim[s], log.get(s)

        # 1. Trimming cannot create reads.
        if int(t["n_reads"]) > int(r["n_reads"]):
            failures.append(
                f"[{s}] trimmed read count {t['n_reads']} exceeds raw "
                f"{r['n_reads']} - trimming cannot create reads"
            )

        # 2. Trimming cannot create bases.
        if int(t["total_bases"]) > int(r["total_bases"]):
            failures.append(
                f"[{s}] trimmed bases {t['total_bases']} exceeds raw "
                f"{r['total_bases']} - trimming cannot create bases"
            )

        # 3. Direction check: if raw and trimmed channels were swapped, this fires.
        if float(t["mean_quality"]) < float(r["mean_quality"]):
            failures.append(
                f"[{s}] mean quality fell after trimming "
                f"({r['mean_quality']} -> {t['mean_quality']}) - "
                f"inputs may be swapped"
            )

        # 4. %Q30 must not decrease.
        if float(t["pct_q30"]) < float(r["pct_q30"]):
            failures.append(
                f"[{s}] pct_q30 fell after trimming "
                f"({r['pct_q30']} -> {t['pct_q30']})"
            )

        # 5. Absolute floor on trimmed output.
        if float(t["pct_q30"]) < min_q30_after:
            failures.append(
                f"[{s}] trimmed pct_q30 {t['pct_q30']} below required floor "
                f"{min_q30_after}"
            )

        # 6. Trim log must reconcile with QC read counts.
        if lg is not None:
            expected = int(lg["reads_in"]) - int(lg["reads_dropped_too_short"])
            if expected != int(lg["reads_out"]):
                failures.append(
                    f"[{s}] trim log inconsistent: reads_in - dropped = "
                    f"{expected} but reads_out = {lg['reads_out']}"
                )
            if int(lg["reads_out"]) != int(t["n_reads"]):
                failures.append(
                    f"[{s}] trim log reads_out {lg['reads_out']} != trimmed "
                    f"QC n_reads {t['n_reads']}"
                )

        # 7. Everything filtered away.
        if int(t["n_reads"]) == 0:
            failures.append(f"[{s}] no reads survived trimming")

    return failures


def main(argv=None):
    ap = argparse.ArgumentParser(description="Validate QC invariants")
    ap.add_argument("--raw-qc", required=True)
    ap.add_argument("--trimmed-qc", required=True)
    ap.add_argument("--trim-log", required=True)
    ap.add_argument("--min-q30-after", type=float, default=70.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    failures = check(
        read_tsv(args.raw_qc),
        read_tsv(args.trimmed_qc),
        read_tsv(args.trim_log),
        args.min_q30_after,
    )

    with open(args.out, "w") as fh:
        if failures:
            fh.write("VALIDATION FAILED\n\n")
            for f in failures:
                fh.write(f"  FAIL  {f}\n")
        else:
            fh.write("VALIDATION PASSED\n\n")
            fh.write("  All QC invariants hold:\n")
            fh.write("    - trimming did not create reads or bases\n")
            fh.write("    - mean quality and %Q30 did not decrease\n")
            fh.write("    - trimmed %Q30 above required floor\n")
            fh.write("    - trim logs reconcile with QC read counts\n")

    if failures:
        print("VALIDATION FAILED:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1

    print("VALIDATION PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
