"""Unit tests for the validation gate — the control that can fail the run."""
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

import validate_qc  # noqa: E402


def row(sample, n_reads=100, total_bases=10000, mean_quality=35.0, pct_q30=90.0):
    return {"sample": sample, "n_reads": str(n_reads),
            "total_bases": str(total_bases),
            "mean_quality": str(mean_quality), "pct_q30": str(pct_q30)}


def log_row(sample, reads_in=100, reads_out=100, dropped=0):
    return {"sample": sample, "reads_in": str(reads_in),
            "reads_out": str(reads_out),
            "reads_dropped_too_short": str(dropped)}


def test_clean_run_passes():
    raw = [row("A", 100, 10000, 30.0, 80.0)]
    trim = [row("A", 95, 8000, 35.0, 92.0)]
    assert validate_qc.check(raw, trim, [log_row("A", 100, 95, 5)]) == []


def test_detects_read_count_increase():
    raw = [row("A", 100)]
    trim = [row("A", 120)]
    fails = validate_qc.check(raw, trim, [log_row("A", 100, 120, 0)])
    assert any("cannot create reads" in f for f in fails)


def test_detects_base_count_increase():
    raw = [row("A", 100, total_bases=10000)]
    trim = [row("A", 100, total_bases=12000)]
    fails = validate_qc.check(raw, trim, [log_row("A")])
    assert any("cannot create bases" in f for f in fails)


def test_detects_swapped_inputs():
    """Passing trimmed as raw must be caught: quality would appear to fall."""
    good = row("A", 95, 8000, 35.0, 92.0)
    bad = row("A", 100, 10000, 30.0, 80.0)
    fails = validate_qc.check([good], [bad], [log_row("A", 95, 100, 0)])
    assert any("inputs may be swapped" in f for f in fails)


def test_detects_q30_regression():
    raw = [row("A", pct_q30=90.0)]
    trim = [row("A", pct_q30=85.0)]
    fails = validate_qc.check(raw, trim, [log_row("A")])
    assert any("pct_q30 fell" in f for f in fails)


def test_enforces_absolute_q30_floor():
    raw = [row("A", pct_q30=40.0, mean_quality=20.0)]
    trim = [row("A", pct_q30=50.0, mean_quality=25.0)]
    fails = validate_qc.check(raw, trim, [log_row("A")], min_q30_after=70.0)
    assert any("below required floor" in f for f in fails)


def test_detects_inconsistent_trim_log():
    raw = [row("A", 100)]
    trim = [row("A", 90)]
    # reads_in - dropped = 95, but reads_out claims 90
    fails = validate_qc.check(raw, trim, [log_row("A", 100, 90, 5)])
    assert any("trim log inconsistent" in f for f in fails)


def test_detects_log_qc_read_count_disagreement():
    raw = [row("A", 100)]
    trim = [row("A", 90)]
    fails = validate_qc.check(raw, trim, [log_row("A", 100, 95, 5)])
    assert any("!= trimmed" in f for f in fails)


def test_detects_total_read_loss():
    raw = [row("A", 100)]
    trim = [row("A", 0, total_bases=0, mean_quality=40.0, pct_q30=100.0)]
    fails = validate_qc.check(raw, trim, [log_row("A", 100, 0, 100)])
    assert any("no reads survived" in f for f in fails)


def test_detects_missing_sample():
    raw = [row("A"), row("B")]
    trim = [row("A")]
    fails = validate_qc.check(raw, trim, [log_row("A")])
    assert any("sample sets differ" in f for f in fails)


def test_cli_exit_codes(tmp_path):
    def write(path, rows):
        cols = list(rows[0].keys())
        Path(path).write_text(
            "\t".join(cols) + "\n"
            + "\n".join("\t".join(r[c] for c in cols) for r in rows) + "\n")

    raw_p, trim_p, log_p = (tmp_path / n for n in ("raw.tsv", "trim.tsv", "log.tsv"))
    write(raw_p, [row("A", 100, 10000, 30.0, 80.0)])
    write(trim_p, [row("A", 95, 8000, 35.0, 92.0)])
    write(log_p, [log_row("A", 100, 95, 5)])
    ok = validate_qc.main(["--raw-qc", str(raw_p), "--trimmed-qc", str(trim_p),
                           "--trim-log", str(log_p), "--out", str(tmp_path / "v.txt")])
    assert ok == 0
    assert "PASSED" in (tmp_path / "v.txt").read_text()

    write(trim_p, [row("A", 120, 20000, 20.0, 50.0)])
    bad = validate_qc.main(["--raw-qc", str(raw_p), "--trimmed-qc", str(trim_p),
                            "--trim-log", str(log_p), "--out", str(tmp_path / "v2.txt")])
    assert bad == 1
    assert "FAILED" in (tmp_path / "v2.txt").read_text()
