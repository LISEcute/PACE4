"""Run complete, independent PACE4 jobs in isolated child processes.

Example:
  python tools/pace_parallel.py --exe _install_MSVC/PACE4.exe \
      --manifest jobs.json --out results/run1 --jobs 4

Manifest: {"version": 1, "jobs": [{"id": "replica-1", "input": "Ni60.pace", "seed": 436}]}
Input paths are relative to the manifest. Every job needs a unique ID and explicit
seed. Inputs are copied unchanged; a calculation is never divided into chunks.
Reproducible scientific files are .cs4, .particles, and the binary traceback;
HTML/RTF also include run-specific dates/paths. Requires Python 3.9 or newer.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone

SCHEMA_VERSION = 1
MAX_CASCADES = 1_000_000
RNG_MODULUS = 2_147_483_647
RNG_MASK = 123_459_876
RESERVED_NAMES = {'CON', 'PRN', 'AUX', 'NUL'} | {f'{x}{i}' for x in ('COM', 'LPT') for i in range(1, 10)}


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def effective_seed(seed):
    state = seed ^ RNG_MASK
    return 1 if state in (0, RNG_MODULUS) else state


def read_manifest(path):
    """Validate and snapshot all inputs before creating outputs or starting children."""
    path = Path(path).resolve()
    raw = path.read_bytes()
    document = json.loads(raw.decode('utf-8-sig'))
    if not isinstance(document, dict) or set(document) != {'version', 'jobs'}:
        raise ValueError('Manifest must contain exactly version and jobs; splitting/NCASC overrides are unsupported.')
    if type(document['version']) is not int or document['version'] != SCHEMA_VERSION:
        raise ValueError('Manifest version must be 1.')
    if not isinstance(document['jobs'], list) or not document['jobs']:
        raise ValueError('Manifest jobs must be a nonempty list.')
    jobs, ids, states = [], set(), set()
    for item in document['jobs']:
        if not isinstance(item, dict) or set(item) != {'id', 'input', 'seed'}:
            raise ValueError('Each job must contain exactly id, input and seed; splitting/overrides are unsupported.')
        name, seed = item['id'], item['seed']
        if not isinstance(name, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', name):
            raise ValueError('Job IDs must be 1-64 ASCII letters, digits, underscores or hyphens, starting with a letter/digit.')
        if name.upper() in RESERVED_NAMES or name.casefold() in ids:
            raise ValueError(f'Duplicate or reserved job ID: {name}')
        if type(seed) is not int or not 1 <= seed < RNG_MODULUS:
            raise ValueError(f'{name}: seed must be an integer from 1 to 2147483646.')
        state = effective_seed(seed)
        if state in states:
            raise ValueError(f'{name}: seed duplicates another effective RNG stream, including degenerate-seed aliases.')
        if not isinstance(item['input'], str) or not item['input']:
            raise ValueError(f'{name}: input must be a nonempty path.')
        source = (path.parent / item['input']).resolve()
        data = source.read_bytes()
        lines = data.decode('utf-8-sig').splitlines()
        fields = lines[0].split() if lines else []
        if len(fields) not in (10, 11) or not re.fullmatch(r'[0-9]+', fields[0]):
            raise ValueError(f'{name}: input must start with a PACE header containing integer NCASC and 10 or 11 fields.')
        ncasc = int(fields[0])
        if not 1 <= ncasc <= MAX_CASCADES:
            raise ValueError(f'{name}: NCASC must be between 1 and {MAX_CASCADES}.')
        ids.add(name.casefold())
        states.add(state)
        jobs.append(dict(id=name, seed=seed, effective_rng_initial_state=state,
                         input_source=str(source), input_sha256=sha256(data), input_bytes=data,
                         ncasc=ncasc, effective_discriminator=1.0 / max(500.0, 3.0 * math.sqrt(ncasc) + 10.0)))
    # Scheduling order cannot affect identity, seed assignment, or output paths.
    jobs.sort(key=lambda job: job['id'])
    return path, raw, jobs


def save_json(path, document):
    temporary = path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(document, indent=2) + '\n', encoding='utf-8')
    temporary.replace(path)


def output_files(folder):
    paths = [p for p in folder.glob('input.*') if p.suffix != '.pace']
    paths += [folder / 'stdout.log', folder / 'stderr.log']
    traces = folder / 'runtime' / 'files'
    if traces.exists():
        paths += list(traces.rglob('*'))
    return {p.relative_to(folder).as_posix(): dict(sha256=sha256(p.read_bytes()), bytes=p.stat().st_size)
            for p in sorted(paths) if p.is_file()}


def child_environment(executable, seed):
    env = dict(os.environ)
    # Development instrumentation must not write to a shared inherited filename.
    for key in list(env):
        if key.startswith('PACE_BENCHMARK_'):
            del env[key]
    env['PACE4_SEED'] = str(seed)
    # Deployed Windows packages often ship qwindows only; silent mode does not show a window.
    offscreen = executable.parent / 'platforms' / 'qoffscreen.dll'
    env['QT_QPA_PLATFORM'] = 'windows' if os.name == 'nt' and not offscreen.is_file() else 'offscreen'
    env['PATH'] = str(executable.parent) + os.pathsep + env.get('PATH', '')
    plugins = [str(executable.parent), str(executable.parent / 'plugins')]
    if env.get('QT_PLUGIN_PATH'):
        plugins.append(env['QT_PLUGIN_PATH'])
    env['QT_PLUGIN_PATH'] = os.pathsep.join(plugins)
    if (executable.parent / 'platforms').is_dir():
        env['QT_QPA_PLATFORM_PLUGIN_PATH'] = str(executable.parent / 'platforms')
    return env


def stop_children(active):
    """Terminate and reap all direct PACE children, escalating after a short grace period."""
    for running in active.values():
        process = running['process']
        if process.poll() is None:
            try:
                if os.name == 'nt':
                    process.terminate()
                else:
                    os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
    deadline = time.monotonic() + 3.0
    for running in active.values():
        process = running['process']
        try:
            process.wait(timeout=max(0.0, deadline - time.monotonic()))
        except subprocess.TimeoutExpired:
            if os.name == 'nt':
                process.kill()
            else:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            process.wait()


def run_suite(executable, manifest, out, workers=1, timeout=None, data_dir=None):
    """Run a manifest; return (exit_code, result). No worker modifies another job's files."""
    started = time.monotonic()
    if type(workers) is not int or workers < 1:
        raise ValueError('--jobs must be a positive integer.')
    if timeout is not None and (not math.isfinite(timeout) or timeout <= 0):
        raise ValueError('--timeout must be a positive finite number of seconds.')
    executable = Path(executable).resolve()
    if not executable.is_file():
        raise ValueError(f'Executable not found: {executable}')
    binary = executable.read_bytes()
    database = (Path(data_dir).resolve() if data_dir else executable.parent / 'lisecfg') / 'AME_DB.sqlite'
    if not database.is_file():
        raise ValueError(f'Mass database not found: {database}; supply --data-dir with the intended lisecfg directory.')
    database_bytes = database.read_bytes()
    manifest, manifest_bytes, jobs = read_manifest(manifest)
    out = Path(out).resolve()
    # Atomic creation refuses existing results, including empty directories.
    out.mkdir(parents=True, exist_ok=False)
    result = dict(version=SCHEMA_VERSION, status='preparing', started_utc=utc_now(),
                  executable=str(executable), executable_sha256=sha256(binary),
                  database=str(database), database_sha256=sha256(database_bytes),
                  manifest=str(manifest), manifest_sha256=sha256(manifest_bytes),
                  workers=workers, timeout_seconds=timeout,
                  seed_policy='Explicit per-job seed; unchanged across worker counts and manifest order.',
                  cutoff_policy='Original NCASC preserved byte-for-byte; no cascade splitting.',
                  runtime_policy='Private executable, database, WorkInROOT ini and traceback directory per job; shared read-only DLL/plugin locations.',
                  setup_seconds=0.0, execution_seconds=0.0, total_seconds=0.0, jobs=[])
    for job in jobs:
        result['jobs'].append({key: value for key, value in job.items() if key != 'input_bytes'})
        result['jobs'][-1].update(status='pending', seconds=0.0, outputs={})
    results_path = out / 'results.json'
    save_json(results_path, result)
    active, next_job, interrupted = {}, 0, []
    execution_started = None
    old_handlers = {}
    old_error_mode = None

    def on_signal(number, frame):
        interrupted.append(number)

    def finish(index, status, error=None):
        running = active.pop(index)
        entry = result['jobs'][index]
        entry.update(status=status, returncode=running['process'].returncode,
                     seconds=time.monotonic()-running['started'], finished_utc=utc_now())
        for stream in running['streams']:
            stream.close()
        if error:
            entry['error'] = error
        entry['outputs'] = output_files(Path(entry['folder']))

    try:
        if os.name == 'nt':
            import ctypes
            # Inherited by children: missing DLLs/crashes must fail, not wait on modal OS dialogs.
            old_error_mode = ctypes.windll.kernel32.SetErrorMode(0x8003)
        # Register before preparation so Ctrl+C never leaves a running child orphaned.
        for number in (signal.SIGINT, signal.SIGTERM):
            old_handlers[number] = signal.signal(number, on_signal)
        (out / 'manifest.json').write_bytes(manifest_bytes)
        for job, entry in zip(jobs, result['jobs']):
            if interrupted:
                break
            folder = out / 'jobs' / job['id']
            runtime = folder / 'runtime'
            (runtime / 'lisecfg').mkdir(parents=True)
            (runtime / 'files').mkdir()
            local_exe = runtime / executable.name
            local_exe.write_bytes(binary)
            if os.name != 'nt':
                local_exe.chmod(executable.stat().st_mode)
            (runtime / 'lisecfg' / 'AME_DB.sqlite').write_bytes(database_bytes)
            (runtime / 'lisepp.ini').write_text('[Version]\nWorkInROOT=1\n[font]\nsize=9\nscaling=0\n', encoding='ascii')
            (folder / 'input.pace').write_bytes(job['input_bytes'])
            entry.update(folder=str(folder), runtime_executable=str(local_exe))
        result['setup_seconds'] = time.monotonic()-started
        execution_started = time.monotonic()
        result['status'] = 'running'
        save_json(results_path, result)
        failed = False
        while (next_job < len(jobs) or active) and not interrupted and not failed:
            while next_job < len(jobs) and len(active) < workers and not interrupted:
                index = next_job
                next_job += 1
                entry = result['jobs'][index]
                folder = Path(entry['folder'])
                streams = [(folder / name).open('wb') for name in ('stdout.log', 'stderr.log')]
                command = [entry['runtime_executable'], 'input.pace', '-s']
                options = dict(cwd=folder, env=child_environment(executable, entry['seed']),
                               stdin=subprocess.DEVNULL, stdout=streams[0], stderr=streams[1])
                if os.name == 'nt':
                    options['creationflags'] = subprocess.CREATE_NO_WINDOW | subprocess.CREATE_NEW_PROCESS_GROUP
                else:
                    options['start_new_session'] = True
                start = time.monotonic()
                try:
                    process = subprocess.Popen(command, **options)
                except OSError as error:
                    for stream in streams:
                        stream.close()
                    entry.update(status='failed', error=str(error), command=command, finished_utc=utc_now())
                    failed = True
                    break
                entry.update(status='running', pid=process.pid, started_utc=utc_now(), command=command)
                active[index] = dict(process=process, started=start, streams=streams)
                save_json(results_path, result)
            for index in list(active):
                running = active[index]
                code = running['process'].poll()
                if code is not None:
                    folder = Path(result['jobs'][index]['folder'])
                    missing = [name for name in ('input.cs4', 'input.html', 'input.rtf')
                               if not (folder / name).is_file() or (folder / name).stat().st_size == 0]
                    if code != 0 or missing:
                        finish(index, 'failed', f'Process exit {code}; missing/empty outputs: {missing}')
                        failed = True
                    else:
                        finish(index, 'succeeded')
                    save_json(results_path, result)
                elif timeout is not None and time.monotonic()-running['started'] >= timeout:
                    result['jobs'][index]['error'] = f'Job exceeded timeout of {timeout} seconds.'
                    result['jobs'][index]['status'] = 'timed_out'
                    failed = True
            if not failed and not interrupted and active:
                time.sleep(0.02)
        if active:
            stop_children(active)
            for index in list(active):
                prior = result['jobs'][index]
                status = 'timed_out' if prior['status'] == 'timed_out' else 'canceled'
                finish(index, status, prior.get('error', 'Stopped after interruption or another job failure.'))
        for entry in result['jobs']:
            if entry['status'] == 'pending':
                entry.update(status='canceled', error='Not started after interruption or another job failure.')
        result['status'] = 'canceled' if interrupted else ('failed' if failed else 'succeeded')
        if interrupted:
            result['signal'] = interrupted[0]
    except BaseException as error:
        stop_children(active)
        for index in list(active):
            finish(index, 'canceled', 'Runner stopped after an internal error.')
        for entry in result['jobs']:
            if entry['status'] == 'pending':
                entry.update(status='canceled', error='Runner stopped before job launch.')
        result.update(status='canceled' if isinstance(error, KeyboardInterrupt) else 'failed', error=str(error))
    finally:
        if old_error_mode is not None:
            ctypes.windll.kernel32.SetErrorMode(old_error_mode)
        for number, handler in old_handlers.items():
            signal.signal(number, handler)
        result['finished_utc'] = utc_now()
        result['total_seconds'] = time.monotonic()-started
        if execution_started is not None:
            result['execution_seconds'] = time.monotonic()-execution_started
        else:
            result['setup_seconds'] = result['total_seconds']
        save_json(results_path, result)
    code = 0 if result['status'] == 'succeeded' else (130 if result['status'] == 'canceled' else 1)
    return code, result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--exe', type=Path, required=True, help='PACE4 executable with deployed runtime DLLs/plugins.')
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True, help='New output directory; existing directories are refused.')
    parser.add_argument('--jobs', type=int, default=1, help='Maximum simultaneous independent processes (default: 1).')
    parser.add_argument('--timeout', type=float, help='Per-job execution timeout in seconds; default: unlimited.')
    parser.add_argument('--data-dir', type=Path, help='Directory containing AME_DB.sqlite; default: executable/lisecfg.')
    args = parser.parse_args()
    try:
        code, result = run_suite(args.exe, args.manifest, args.out, args.jobs, args.timeout, args.data_dir)
    except (OSError, ValueError) as error:
        print(f'PACE4 parallel runner: {error}', file=sys.stderr)
        return 2
    print(f"{result['status']}: {len(result['jobs'])} jobs, {result['total_seconds']:.3f}s total "
          f"({result['setup_seconds']:.3f}s setup); {args.out.resolve() / 'results.json'}")
    for job in result['jobs']:
        if job['status'] != 'succeeded':
            print(f"  {job['id']}: {job['status']}: {job.get('error', '')}", file=sys.stderr)
    return code


if __name__ == '__main__':
    sys.exit(main())
