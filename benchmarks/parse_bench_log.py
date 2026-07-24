#!/usr/bin/env python3
"""
parse_bench_log.py <pasa_run.log>

Parses Pipeliner's verbose "* [<localtime>] Running CMD: <cmdstr>" lines from
a Launch_PASA_pipeline.pl run log, computes the wall-clock gap between
consecutive commands, and aggregates by canonical sub-script name (basename
of the first token of the command line). This directly measures each PASA
sub-step's real duration for a given run -- used as the before/after
benchmark for the assign_clusters_by_*.dbi / subcluster_builder.dbi /
TransDecoder threading fixes.
"""
import re
import statistics
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

LINE_RE = re.compile(r'^\* \[(.*?)\] Running CMD: (.*)$')
# Perl localtime() stringifies as e.g. "Wed Jul 23 14:32:10 2026"
TS_FMT = '%a %b %d %H:%M:%S %Y'


def canonical_name(cmdstr):
    tok = cmdstr.strip().split()[0]
    return Path(tok).name


def main():
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <pasa_run.log>", file=sys.stderr)
        sys.exit(1)

    log_path = Path(sys.argv[1])
    events = []
    with open(log_path, errors='replace') as fh:
        for line in fh:
            m = LINE_RE.match(line.rstrip('\n'))
            if not m:
                continue
            raw_ts, cmdstr = m.group(1), m.group(2)
            # normalize double-space single-digit day ("Jul  3") -> strptime tolerant
            ts = datetime.strptime(' '.join(raw_ts.split()), TS_FMT)
            events.append((ts, canonical_name(cmdstr), cmdstr))

    if len(events) < 2:
        print(f"Only found {len(events)} 'Running CMD' lines in {log_path}; nothing to report.",
              file=sys.stderr)
        sys.exit(1)

    durations = defaultdict(list)
    for i in range(len(events) - 1):
        ts0, name0, _ = events[i]
        ts1, _, _ = events[i + 1]
        dur = (ts1 - ts0).total_seconds()
        if dur < 0:
            continue
        durations[name0].append(dur)

    total_wall = (events[-1][0] - events[0][0]).total_seconds()

    rows = []
    for name, durs in durations.items():
        rows.append((name, len(durs), sum(durs), statistics.median(durs), max(durs)))
    rows.sort(key=lambda r: -r[2])

    print(f"\n{'Script':<45} {'N':<4} {'Total (s)':<12} {'Median (s)':<12} {'Max (s)':<10} {'% of measured':<10}")
    print('-' * 100)
    grand_total = sum(r[2] for r in rows) or 1
    for name, n, total, median, mx in rows:
        print(f"{name:<45} {n:<4} {total:<12.1f} {median:<12.1f} {mx:<10.1f} {100*total/grand_total:<10.1f}")
    print('-' * 100)
    print(f"Total measured (sum of gaps): {grand_total:.1f}s | "
          f"First-to-last-CMD wall time: {total_wall:.1f}s\n")


if __name__ == '__main__':
    main()
