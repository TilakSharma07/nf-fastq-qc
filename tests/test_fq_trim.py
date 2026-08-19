"""Unit tests for the FASTQ trimmer."""
import gzip
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

import fq_qc  # noqa: E402
import fq_trim  # noqa: E402

ADAPTER = "AGATCGGAAGAGC"


def write_fastq(path, records):
    Path(path).write_text("".join(f"@{h}\n{s}\n+\n{q}\n" for h, s, q in records))
    return path


def read_fastq_gz(path):
    with gzip.open(path, "rt") as fh:
        return list(fq_qc.iter_fastq(fh))


# ---------------- adapter matching ----------------

def test_find_adapter_exact_match():
    seq = "ACGTACGT" + ADAPTER
    assert fq_trim.find_adapter(seq, ADAPTER) == 8


def test_find_adapter_tolerates_one_mismatch():
    mutated = "T" + ADAPTER[1:]           # 1/13 mismatches = 7.7% < 15%
    assert fq_trim.find_adapter("AAAA" + mutated, ADAPTER) == 4


def test_find_adapter_rejects_excess_mismatch():
    junk = "TTTTTTTTTTTTT"
    assert fq_trim.find_adapter("AAAA" + junk, ADAPTER) == -1


def test_find_adapter_absent_returns_minus_one():
    assert fq_trim.find_adapter("ACGTACGTACGTACGT", ADAPTER) == -1


def test_find_adapter_partial_at_read_end():
    # only the first 6 adapter bases are present, at the 3' end
    seq = "ACGTACGTAC" + ADAPTER[:6]
    assert fq_trim.find_adapter(seq, ADAPTER) == 10


def test_find_adapter_respects_min_overlap():
    seq = "ACGTACGTAC" + ADAPTER[:2]      # 2 bp overlap < min_overlap 3
    assert fq_trim.find_adapter(seq, ADAPTER, min_overlap=3) == -1


# ---------------- quality trimming ----------------

def test_trim_trailing_quality_removes_bad_tail():
    seq, qual = "ACGTACGT", "IIII####"    # Phred 40 x4 then 2 x4
    s, q = fq_trim.trim_trailing_quality(seq, qual, 20)
    assert s == "ACGT"
    assert q == "IIII"
    assert len(s) == len(q)


def test_trim_trailing_quality_keeps_good_read_intact():
    seq, qual = "ACGTACGT", "IIIIIIII"
    s, q = fq_trim.trim_trailing_quality(seq, qual, 20)
    assert s == seq and q == qual


def test_trim_trailing_quality_tolerates_isolated_dip():
    # a single Phred-2 base inside a Phred-40 tail must not truncate the read
    qual = "IIII" + "#" + "IIII"
    s, _ = fq_trim.trim_trailing_quality("A" * 9, qual, 20)
    assert len(s) == 9


def test_trim_trailing_quality_can_empty_a_read():
    s, q = fq_trim.trim_trailing_quality("ACGT", "####", 20)
    assert s == "" and q == ""


def test_sequence_and_quality_stay_in_lockstep():
    s, q = fq_trim.trim_trailing_quality("ACGTACGTAC", "IIIII#####", 20)
    assert len(s) == len(q)


# ---------------- end-to-end ----------------

def test_min_length_filter_drops_short_reads(tmp_path):
    fq = write_fastq(tmp_path / "in.fastq", [
        ("keep", "A" * 60, "I" * 60),
        ("drop", "A" * 60, "I" * 10 + "#" * 50),   # trims to 10 bp
    ])
    out, log = tmp_path / "o.fastq.gz", tmp_path / "o.log.tsv"
    rc = fq_trim.main(["--fastq", str(fq), "--sample", "S",
                       "--out-fastq", str(out), "--out-log", str(log),
                       "--quality-cutoff", "20", "--min-length", "25"])
    assert rc == 0
    kept = read_fastq_gz(out)
    assert len(kept) == 1
    assert kept[0][0] == "@keep"
    row = dict(zip(*[l.split("\t") for l in log.read_text().strip().split("\n")]))
    assert row["reads_in"] == "2"
    assert row["reads_out"] == "1"
    assert row["reads_dropped_too_short"] == "1"


def test_trimming_never_increases_bases(tmp_path):
    """The core invariant the pipeline's validation gate also enforces."""
    fq = write_fastq(tmp_path / "in.fastq",
                     [(f"r{i}", "ACGT" * 25, "I" * 60 + "#" * 40) for i in range(20)])
    out, log = tmp_path / "o.fastq.gz", tmp_path / "o.log.tsv"
    fq_trim.main(["--fastq", str(fq), "--sample", "S", "--out-fastq", str(out),
                  "--out-log", str(log)])
    before = fq_qc.compute_metrics(fq, "S")
    after = fq_qc.compute_metrics(out, "S")
    assert after["total_bases"] <= before["total_bases"]
    assert after["n_reads"] <= before["n_reads"]


def test_trimming_improves_mean_quality(tmp_path):
    fq = write_fastq(tmp_path / "in.fastq",
                     [(f"r{i}", "ACGT" * 25, "I" * 60 + "#" * 40) for i in range(20)])
    out, log = tmp_path / "o.fastq.gz", tmp_path / "o.log.tsv"
    fq_trim.main(["--fastq", str(fq), "--sample", "S", "--out-fastq", str(out),
                  "--out-log", str(log)])
    before = fq_qc.compute_metrics(fq, "S")
    after = fq_qc.compute_metrics(out, "S")
    assert after["mean_quality"] > before["mean_quality"]
    assert after["pct_q30"] > before["pct_q30"]


def test_output_is_deterministic(tmp_path):
    fq = write_fastq(tmp_path / "in.fastq",
                     [(f"r{i}", "ACGT" * 25, "I" * 55 + "#" * 45) for i in range(30)])
    hashes = []
    for k in range(2):
        out = tmp_path / f"o{k}.fastq.gz"
        fq_trim.main(["--fastq", str(fq), "--sample", "S", "--out-fastq", str(out),
                      "--out-log", str(tmp_path / f"l{k}.tsv")])
        hashes.append(read_fastq_gz(out))
    assert hashes[0] == hashes[1]


def test_trimmed_output_is_valid_fastq(tmp_path):
    fq = write_fastq(tmp_path / "in.fastq",
                     [("r1", "ACGT" * 25 , "I" * 70 + "#" * 30)])
    out = tmp_path / "o.fastq.gz"
    fq_trim.main(["--fastq", str(fq), "--sample", "S", "--out-fastq", str(out),
                  "--out-log", str(tmp_path / "l.tsv")])
    # iter_fastq raises on any structural problem, so parsing is the assertion
    recs = read_fastq_gz(out)
    assert all(len(s) == len(q) for _h, s, q in recs)
