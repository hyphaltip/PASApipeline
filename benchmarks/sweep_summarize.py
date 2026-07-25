#!/usr/bin/env python3
"""Summarize a sweep_threads.sh results TSV into a scaling table.

Reports, per script, wall time and speedup relative to the lowest thread count
measured, plus parallel efficiency and peak RSS growth. Flags any thread count
whose output checksum diverges from the serial reference -- a thread count that
changes results is a correctness bug, not a speedup.

Usage: sweep_summarize.py sweep_results.tsv
"""
import sys
import csv
from collections import defaultdict


def main(path):
    rows = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows[r["script"]].append(r)

    for script, rs in rows.items():
        rs.sort(key=lambda r: int(r["threads"]))
        base = rs[0]
        base_t = float(base["wall_s"])
        base_n = int(base["threads"])
        base_ck = base["output_checksum"]

        print()
        print(f"=== {script} (baseline: -T {base_n}) ===")
        print(f"{'threads':>8} {'wall_s':>10} {'speedup':>9} {'efficiency':>11} "
              f"{'peak_RSS_MB':>12} {'output':>9}")
        print("-" * 66)
        for r in rs:
            n = int(r["threads"])
            t = float(r["wall_s"])
            speedup = base_t / t if t else float("nan")
            # efficiency vs the extra parallelism actually added over baseline
            eff = speedup / (n / base_n) if n else float("nan")
            rss = float(r["max_rss_kb"] or 0) / 1024
            if r["exit"] != "0":
                status = f"FAIL({r['exit']})"
            elif r["output_checksum"] != base_ck:
                status = "DIVERGED"
            else:
                status = "match"
            print(f"{n:>8} {t:>10.1f} {speedup:>8.2f}x {eff:>10.0%} "
                  f"{rss:>12.1f} {status:>9}")

        best = min(rs, key=lambda r: float(r["wall_s"]))
        print(f"\n  fastest: -T {best['threads']} at {float(best['wall_s']):.1f}s "
              f"({base_t / float(best['wall_s']):.2f}x vs -T {base_n})")
        diverged = [r["threads"] for r in rs if r["output_checksum"] != base_ck]
        if diverged:
            print(f"  WARNING: output differs from -T {base_n} at threads: "
                  f"{', '.join(diverged)}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
