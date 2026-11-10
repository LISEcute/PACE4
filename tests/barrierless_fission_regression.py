"""End-to-end checks for rotational instability, sampling, and yield accounting.

Run with --exe, --baseline (4.38.3), and a new --out directory. Both executables
need AME_DB.sqlite beside them in lisecfg and their Qt/compiler runtimes on PATH.
The 48Ca+124Sn case reproduces the reported seed and 10,000-event calculation.
"""
import argparse
import collections
import html
import json
from pathlib import Path
import re
import struct
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from pace_parallel import run_suite


def residual_counts(report):
    table = report.split('1. Yields of residual nuclei', 1)[1].split('</table>', 1)[0]
    counts = {}
    for row in re.findall(r'<tr[^>]*>(.*?)</tr>', table, re.S):
        cells = [html.unescape(re.sub('<[^>]+>', '', cell)).strip()
                 for cell in re.findall(r'<td[^>]*>(.*?)</td>', row, re.S)]
        if len(cells) == 7 and all(re.fullmatch(r'\d+', cell) for cell in cells[:3]):
            counts[f'{cells[0]},{cells[1]}'] = int(cells[4])
        elif cells and 'Total fission' in cells[0]:
            counts['fission'] = int(cells[2])
    return counts


def fusion(events, fyrst=0, distributions=1):
    return (f'{events} {distributions} 1 {fyrst} 0 1 0 10 0 0 0\n'
            '20 48 50 124 0 0 0\n300 0 0 2 0 0 0 3\n\n\n\n\n')


def fixed(z, a, energy, spin, fyrst=0):
    return f'64 1 2 {fyrst} 0 1 0 10 0 0 0\n{z} {a} {energy} 0 {spin} 0\n\n\n\n'


def read_job(folder):
    report = (folder / 'input.html').read_text(encoding='utf-8')
    counts = residual_counts(report)
    match = re.search(r'Barrierless fission \(included in Total fission\): (\d+) events', report)
    barrierless = int(match[1]) if match else 0
    events = collections.defaultdict(list)
    trace = folder / 'runtime/files/PACE4evt.dat'
    if trace.exists():
        for record in struct.iter_unpack('<10h5f', trace.read_bytes()):
            event_id = (record[0] << 16) + (record[1] & 65535)
            events[event_id].append(record)
    return report, counts, barrierless, events


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', required=True, type=Path)
    parser.add_argument('--baseline', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    root = args.out.resolve()
    root.mkdir(parents=True, exist_ok=False)
    cases = {
        'exact': (fusion(10000), 708849395),
        'replica': (fusion(1000), 436),
        'cps': (fusion(1000, fyrst=1), 104729),
        'no_trace': (fusion(64, distributions=0), 130363),
        'boundary_below': (fixed(70, 172, 143.076, 102), 155921),
        'boundary_at': (fixed(70, 172, 143.076, 103), 196613),
        'high_spin': (fixed(70, 172, 143.076, 117), 301342),
        'low_spin': (fixed(28, 60, 30.125, 10), 406071),
        'low_spin_odd': (fixed(28, 61, 30.125, 10.5), 510800),
        'bound_gamma_trap': (fixed(70, 172, 1, 10), 615529),
    }
    jobs = []
    for name, (data, seed) in cases.items():
        (root / f'{name}.pace').write_text(data, encoding='ascii')
        jobs.append(dict(id=name, input=f'{name}.pace', seed=seed))
    manifest = root / 'jobs.json'
    manifest.write_text(json.dumps(dict(version=1, jobs=jobs)))
    for label, exe in [('baseline', args.baseline), ('corrected', args.exe)]:
        code, result = run_suite(exe.resolve(), manifest, root / label, 2, 180)
        assert code == 0, (label, result)

    evidence = {}
    for name in cases:
        old_folder = root / 'baseline/jobs' / name
        new_folder = root / 'corrected/jobs' / name
        _, old_counts, _, old_events = read_job(old_folder)
        report, counts, barrierless, events = read_job(new_folder)
        expected = int(cases[name][0].split()[0])
        assert sum(counts.values()) == expected, (name, counts)
        assert barrierless <= counts.get('fission', 0), name
        assert '"FIND" limits error' not in report, name
        assert 'LEVEL DENSITY EXCEEDED' not in report, name
        if name != 'no_trace':
            assert len(events) == len(old_events) == expected, name
            # Compare each sampled event, not just its histogram or mean spin.
            assert {i: r[0][3] for i, r in events.items()} == {
                i: r[0][3] for i, r in old_events.items()}, name
            assert sum(records[-1][2] == 6 for records in events.values()) == counts.get('fission', 0)
            for records in events.values():
                first = records[0]
                if first[8:10] == (70, 102) and first[3]-1 >= 103:
                    assert len(records) == 1 and first[2] == 6, (name, first)
                    assert first[10] == 1 and first[12:15] == (0., 0., 0.), (name, first)
        if name in ('exact', 'replica', 'cps', 'no_trace', 'boundary_at', 'high_spin'):
            assert counts.get('70,102', 0) == 0, (name, counts)
        if name in ('boundary_at', 'high_spin'):
            assert barrierless == counts['fission'] == expected, name
        if name in ('low_spin', 'low_spin_odd', 'bound_gamma_trap'):
            assert barrierless == 0 and counts == old_counts, name
            for filename in ('input.cs4', 'runtime/files/PACE4evt.dat'):
                assert (new_folder / filename).read_bytes() == (old_folder / filename).read_bytes(), (name, filename)
        if name == 'exact':
            assert old_counts.get('70,102') == 2446, 'Baseline must reproduce reported 4.38.3 run'
        evidence[name] = dict(events=expected, before=old_counts, after=counts,
                              barrierless_fission=barrierless)
    (root / 'verification.json').write_text(json.dumps(evidence, indent=2) + '\n')
    print('PASS 10 cases: unchanged sampled spins, conserved events, barrierless endpoints, '
          'and byte-identical low-spin/gamma-trap controls')


if __name__ == '__main__':
    main()
