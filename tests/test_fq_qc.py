"""Unit tests for the FASTQ QC metrics implementation."""
import gzip
import io
import json
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

import fq_qc  # noqa: E402


def write_fastq(path, records, gzipped=False):
    text = "".join(f"@{h}\n{s}\n+\n{q}\n" for h, s, q in records)
    if gzipped:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        Path(path).write_text(text)
    return path


# ---------------- parsing ----------------

def test_iter_fastq_parses_all_records():
    recs = [("r1", "ACGT", "IIII"), ("r2", "TTTT", "JJJJ")]
    fh = io.StringIO("".join(f"@{h}\n{s}\n+\n{q}\n" for h, s, q in recs))
    parsed = list(fq_qc.iter_fastq(fh))
    assert len(parsed) == 2
    assert parsed[0] == ("@r1", "ACGT", "IIII")


def test_iter_fastq_rejects_length_mismatch():
    fh = io.StringIO("@r1\nACGTA\n+\nIIII\n")
    with pytest.raises(ValueError, match="length mismatch"):
        list(fq_qc.iter_fastq(fh))


def test_iter_fastq_rejects_bad_header():
    fh = io.StringIO("r1\nACGT\n+\nIIII\n")
    with pytest.raises(ValueError, match="does not start with"):
        list(fq_qc.iter_fastq(fh))


def test_iter_fastq_rejects_missing_plus():
    fh = io.StringIO("@r1\nACGT\nX\nIIII\n")
    with pytest.raises(ValueError, match="separator"):
        list(fq_qc.iter_fastq(fh))


def test_gzip_and_plain_agree(tmp_path):
    recs = [("r1", "ACGTACGT", "IIIIIIII"), ("r2", "GGCCGGCC", "JJJJJJJJ")]
    plain = write_fastq(tmp_path / "a.fastq", recs, gzipped=False)
    gz = write_fastq(tmp_path / "a.fastq.gz", recs, gzipped=True)
    assert fq_qc.compute_metrics(plain, "s") == fq_qc.compute_metrics(gz, "s")


# ---------------- Phred conversion ----------------

def test_phred_scores_known_values():
    # '!' is Phred 0, 'I' is Phred 40 under the +33 offset.
    assert fq_qc.phred_scores("!") == [0]
    assert fq_qc.phred_scores("I") == [40]
    assert fq_qc.phred_scores("#5I") == [2, 20, 40]


# ---------------- metrics ----------------

def test_gc_percent_exact(tmp_path):
    # 4 G/C out of 8 non-N bases -> 50%
    p = write_fastq(tmp_path / "gc.fastq", [("r1", "GGCCATAT", "IIIIIIII")])
    m = fq_qc.compute_metrics(p, "s")
    assert m["gc_percent"] == pytest.approx(50.0)


def test_gc_percent_excludes_ambiguous(tmp_path):
    # GC=2, AT=2, N=4 -> 50%, N excluded from the denominator
    p = write_fastq(tmp_path / "n.fastq", [("r1", "GCATNNNN", "IIIIIIII")])
    m = fq_qc.compute_metrics(p, "s")
    assert m["gc_percent"] == pytest.approx(50.0)
    assert m["n_bases_ambiguous"] == 4


def test_q20_q30_fractions(tmp_path):
    # Phred 40,40,20,2 -> Q20+: 3/4 = 75%, Q30+: 2/4 = 50%
    p = write_fastq(tmp_path / "q.fastq", [("r1", "ACGT", "II5#")])
    m = fq_qc.compute_metrics(p, "s")
    assert m["pct_q20"] == pytest.approx(75.0)
    assert m["pct_q30"] == pytest.approx(50.0)
    assert m["mean_quality"] == pytest.approx((40 + 40 + 20 + 2) / 4)


def test_read_length_stats(tmp_path):
    p = write_fastq(tmp_path / "len.fastq",
                    [("r1", "AC", "II"), ("r2", "ACGT", "IIII")])
    m = fq_qc.compute_metrics(p, "s")
    assert m["n_reads"] == 2
    assert m["total_bases"] == 6
    assert m["min_read_length"] == 2
    assert m["max_read_length"] == 4
    assert m["mean_read_length"] == pytest.approx(3.0)


def test_per_cycle_quality_handles_ragged_lengths(tmp_path):
    # cycle 1: mean(40, 20) = 30 ; cycle 2: only the longer read contributes
    p = write_fastq(tmp_path / "rag.fastq",
                    [("r1", "A", "I"), ("r2", "AC", "5I")])
    m = fq_qc.compute_metrics(p, "s")
    assert m["per_cycle_mean_quality"][0] == pytest.approx(30.0)
    assert m["per_cycle_mean_quality"][1] == pytest.approx(40.0)
    assert len(m["per_cycle_mean_quality"]) == 2


def test_empty_fastq_raises(tmp_path):
    p = tmp_path / "empty.fastq"
    p.write_text("")
    with pytest.raises(ValueError, match="no FASTQ records"):
        fq_qc.compute_metrics(p, "s")


# ---------------- CLI ----------------

def test_cli_writes_tsv_and_json(tmp_path):
    fq = write_fastq(tmp_path / "cli.fastq", [("r1", "GGCC", "IIII")])
    tsv, js = tmp_path / "o.tsv", tmp_path / "o.json"
    rc = fq_qc.main(["--fastq", str(fq), "--sample", "S1",
                     "--out-tsv", str(tsv), "--out-json", str(js)])
    assert rc == 0
    lines = tsv.read_text().strip().split("\n")
    assert lines[0].split("\t") == fq_qc.TSV_COLUMNS
    assert lines[1].split("\t")[0] == "S1"
    assert json.loads(js.read_text())["sample"] == "S1"


def test_cli_is_executable_as_script(tmp_path):
    fq = write_fastq(tmp_path / "x.fastq", [("r1", "GGCC", "IIII")])
    r = subprocess.run(
        [sys.executable, str(BIN / "fq_qc.py"), "--fastq", str(fq),
         "--sample", "S", "--out-tsv", str(tmp_path / "t.tsv"),
         "--out-json", str(tmp_path / "t.json")],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr


def test_per_cycle_read_count_matches_lengths(tmp_path):
    """Cycle support counts must equal the number of reads reaching that cycle."""
    p = tmp_path / "r.fastq"
    p.write_text(
        "@a\nACGT\n+\nIIII\n"
        "@b\nACG\n+\nIII\n"
        "@c\nAC\n+\nII\n"
    )
    m = fq_qc.compute_metrics(str(p), "s")
    assert m["per_cycle_read_count"] == [3, 3, 2, 1]
    assert len(m["per_cycle_mean_quality"]) == len(m["per_cycle_read_count"])
