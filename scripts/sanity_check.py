#!/usr/bin/env python3
"""
Parses the combined build log written by org.reactome.fi.PipelineRunner (logs/build_*.log)
for the checkpoint metrics the manual procedure treats as review checkpoints (ids
merged per data source, FI counts per source, threshold-vs-count numbers, final
totals), and diffs them against the previous year's saved baseline. Replaces
manually eyeballing the raw log for surprising numbers.

Usage:
    python3 scripts/sanity_check.py [path/to/build.log]
"""
import glob
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOGS_DIR = os.path.join(ROOT, "logs")

# Matches lines like "Total ids from KEGG: 1234" or "Threshold 0.05..." (log4j
# prefixes like "2024-04-16 ... INFO  org.reactome.fi.X  - " are stripped first).
METRIC_RE = re.compile(r'^\s*([A-Za-z][\w /\-]*?):\s*(-?\d+(?:\.\d+)?)\s*(?:\(|$)')
SECTION_RE = re.compile(r'Done data source:\s*(.+)$')
FLAG_THRESHOLD = 0.15


def latest_log():
    logs = sorted(glob.glob(os.path.join(LOGS_DIR, "build_*.log")))
    if not logs:
        sys.exit("No build logs found under logs/ - run scripts/run_pipeline.sh first.")
    return logs[-1]


def strip_log4j_prefix(line):
    # log4j pattern is "%d{ISO8601} [%t] %-5p %c %x - %m%n"; keep just the message.
    if " - " in line:
        return line.split(" - ", 1)[-1]
    return line


def parse_log(path):
    metrics = {}
    section = None
    with open(path, errors="replace") as f:
        for raw_line in f:
            m = SECTION_RE.search(raw_line)
            if m:
                section = m.group(1).strip()
                continue
            m = METRIC_RE.search(strip_log4j_prefix(raw_line))
            if m:
                label, value = m.group(1).strip(), float(m.group(2))
                key = "%s :: %s" % (section, label) if section else label
                metrics[key] = value
    return metrics


def main():
    year = os.environ.get("YEAR")
    log_path = sys.argv[1] if len(sys.argv) > 1 else latest_log()
    metrics = parse_log(log_path)
    print("Parsed %d checkpoint metrics from %s" % (len(metrics), log_path))

    if not year:
        print("YEAR not set in the environment (source .build_env.sh first) - skipping baseline diff.")
        return

    baseline_path = os.path.join(LOGS_DIR, "sanity_baseline_%s.json" % year)
    prev_year = str(int(year) - 1)
    prev_baseline_path = os.path.join(LOGS_DIR, "sanity_baseline_%s.json" % prev_year)

    if os.path.exists(prev_baseline_path):
        with open(prev_baseline_path) as f:
            prev = json.load(f)
        print("\nDeltas vs %s (>=%.0f%% change or a metric that disappeared):" % (prev_year, FLAG_THRESHOLD * 100))
        flagged = False
        for key, value in sorted(metrics.items()):
            if key in prev and prev[key]:
                delta = (value - prev[key]) / abs(prev[key])
                if abs(delta) >= FLAG_THRESHOLD:
                    flagged = True
                    print("  %s: %g -> %g (%+.0f%%)" % (key, prev[key], value, delta * 100))
        for key in prev:
            if key not in metrics:
                flagged = True
                print("  %s: %g -> MISSING" % (key, prev[key]))
        if not flagged:
            print("  (no anomalies)")
    else:
        print("\nNo baseline found for %s - this run's numbers become next year's baseline." % prev_year)

    os.makedirs(LOGS_DIR, exist_ok=True)
    with open(baseline_path, "w") as f:
        json.dump(metrics, f, indent=2, sort_keys=True)
    print("\nSaved baseline to %s" % baseline_path)

    print("\nAlso glance at the compareFilesToPreviousVersion() output further up this log")
    print("(a file-size-level comparison table) for a complementary regression check.")


if __name__ == "__main__":
    main()
