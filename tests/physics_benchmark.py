"""Run reproducible corrected-PACE4 physics benchmark cases and preserve raw reports."""
import argparse
import csv
import hashlib
import html
import json
import os
from pathlib import Path
import re
import subprocess
import time

SEEDS = [436, 104729, 130363, 155921, 196613]
CASES = [
    dict(name="ni60_fixed_j10", mode=2, z=28, a=60, energy=30.125, spin=10),
    dict(name="ni60_equal_j0_10_20", mode=3, z=28, a=60, energy=30.125, spin=20),
    dict(name="ni60_generated_j0_20", mode=5, z=28, a=60, energy=30.125, spin=20),
]


def input_text(case, events, fyrst=0):
    lines = [f"{events} 0 {case['mode']} {fyrst} -1 1 0 10 0 1 0",
             f"{case['z']} {case['a']} {case['energy']} 0 {case['spin']} 0"]
    if case['mode'] == 3:
        weights = [100 / 3 if j in (0, 10, 20) else 0 for j in range(21)]
        lines.append("21 " + " ".join(format(w, '.17g') for w in weights))
    return "\n".join(lines + ["", "", ""]) + "\n"


def residual_counts(report):
    table = report.split('1. Yields of residual nuclei', 1)[1].split('</table>', 1)[0]
    counts = {}
    for row in re.findall(r'<tr[^>]*>(.*?)</tr>', table, re.S):
        cells = [html.unescape(re.sub('<[^>]+>', '', cell)).strip()
                 for cell in re.findall(r'<td[^>]*>(.*?)</td>', row, re.S)]
        if len(cells) == 7 and all(re.fullmatch(r'\d+', cell) for cell in cells[:3]):
            counts[f"{cells[0]},{cells[1]}"] = int(cells[4])
        elif cells and 'Total fission' in cells[0]:
            counts['fission'] = int(cells[2])
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exe', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--events', type=int, default=5000)
    ap.add_argument('--fyrst', type=float, default=0.)
    ap.add_argument('--no-spin-recording', action='store_true', help='Use an uninstrumented production executable')
    ap.add_argument('--replicas', type=int, default=5, choices=range(1, 6))
    args = ap.parse_args()
    exe = args.exe.resolve()
    root = args.out.resolve()
    root.mkdir(parents=True, exist_ok=True)
    suite = dict(executable=str(exe), executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
                 events_per_run=args.events, fyrst=args.fyrst, noshl=1, facla=10, aratio=1, barfac=-1, runs=[])
    for case in CASES:
        for seed in SEEDS[:args.replicas]:
            folder = root / case['name'] / str(seed)
            folder.mkdir(parents=True, exist_ok=True)
            source = folder / 'input.pace'
            source.write_text(input_text(case, args.events, args.fyrst), encoding='ascii')
            env = dict(os.environ, QT_QPA_PLATFORM='offscreen', PACE4_SEED=str(seed),
                       PACE_BENCHMARK_SPINS=str(folder / 'spins.csv'))
            start = time.monotonic()
            result = subprocess.run([str(exe), str(source), '-s'], cwd=folder, env=env,
                                    capture_output=True, timeout=180)
            duration = time.monotonic() - start
            (folder / 'stdout.log').write_bytes(result.stdout)
            (folder / 'stderr.log').write_bytes(result.stderr)
            if result.returncode != 0:
                raise RuntimeError(f"{case['name']} seed {seed}: exit {result.returncode}; inspect {folder}")
            report = source.with_suffix('.html').read_text(encoding='utf-8')
            assert 'Version 4.37.0' in report, 'Benchmark requires PACE4 4.37.0'
            counts = residual_counts(report)
            assert sum(counts.values()) == args.events, (case['name'], seed, counts)
            spins = {}
            if (folder / 'spins.csv').exists():
                with (folder / 'spins.csv').open() as f:
                    spins = {int(row['spin_index'])-1:int(row['count']) for row in csv.DictReader(f)}
            if case['mode'] == 2:
                spins = {case['spin']:args.events}
            elif not args.no_spin_recording:
                assert sum(spins.values()) == args.events, (case['name'], seed, spins)
            cs_header = source.with_suffix('.cs4').read_text().splitlines()[1].split()
            run = dict(case=case, seed=seed, seconds=duration, events=args.events,
                       sigma_input_mb=float(cs_header[4]), residual_counts=counts,
                       initial_spin_counts=spins, folder=str(folder))
            suite['runs'].append(run)
            (root / 'results.json').write_text(json.dumps(suite, indent=2))
            print(f"{case['name']} seed={seed}: {args.events} events, {duration:.2f}s", flush=True)


if __name__ == '__main__':
    main()
