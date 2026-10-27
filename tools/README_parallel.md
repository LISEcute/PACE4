# Independent PACE4 processes

## Run from the PACE4 interface

Choose **Parallel calculations…** from the File menu or toolbar. Python is not required for this mode.

1. Save each calculation as a `.pace` file using the ordinary PACE4 input editor.
2. In the parallel dialog, add the saved files. Duplicate a row to run another replica of the same input.
3. Set each row's seed and choose the number of workers. Seeds stay attached to jobs when the worker count changes.
4. Choose a new results folder, then start the batch. Progress and per-job status appear in the dialog.
5. Select a completed row and use **Open results** to view its HTML report. With no row selected, this opens the batch folder. Use **Cancel batch** to stop active and pending jobs.

Only saved input files are used; unsaved edits in the main window are not included. Each job retains its full cascade count and probability cutoff. Changing the worker count changes concurrency, not the scientific inputs. Different seeds produce different replicas; keep the same input and seed to replay a calculation.

Each process has its own executable copy, database, configuration and traceback directory. The batch records input and output hashes, seeds, cascade counts, cutoffs and status in `results.json`. Keep the PACE4 runtime libraries and plugins unchanged during a batch. No extra Python installation or command-line setup is needed.

## Command-line runner

`pace_parallel.py` runs complete, independent calculations concurrently. It uses Python 3.9 or newer and the standard library. The PACE4 installation must include its matching runtime libraries, platform/SQL plugins, and `lisecfg/AME_DB.sqlite`.

From the repository root, run the included example:

```powershell
python tools/pace_parallel.py --exe _install_MSVC/PACE4.exe --manifest tools/examples/parallel/jobs.json --out test/parallel_example --jobs 4 --timeout 180
```

Use `_install/PACE4.exe` for the MinGW installation. Choose a **new output directory** for every run. The runner refuses an existing directory rather than mixing old and new results. `--jobs 1` is the serial reference, using exactly the same isolation and seed policy as parallel runs. `--timeout` limits each child calculation, not the whole batch.

## Manifest

```json
{
  "version": 1,
  "jobs": [
    {"id": "fixed-1", "input": "fixed.pace", "seed": 436},
    {"id": "fixed-2", "input": "fixed.pace", "seed": 104729}
  ]
}
```

Input paths are relative to the manifest. Each job needs a unique ID and an explicit integer seed from 1 through 2147483646. IDs start with an ASCII letter or digit and contain only ASCII letters, digits, underscores or hyphens; Windows reserved names are rejected. Duplicate effective RNG starting states, including masked-seed aliases, are rejected. Explicit distinct starting states do not establish mathematically nonoverlapping random streams for arbitrary run lengths.

The ID, input and seed identify the calculation. They do not depend on scheduling or the number of workers. The runner snapshots input bytes, the executable and the mass database before starting. Preserve the executable, database, manifest and inputs to reproduce a run.

## Preserve the probability cutoff

The input file is copied byte-for-byte. NCASC is neither overridden nor divided among workers. PACE4 calculates:

```text
DISCPR = 1 / max(500, 3 * sqrt(NCASC) + 10)
```

Every manifest entry remains one complete calculation. The runner does not accept cascade-splitting options. Repeating a complete job with another seed is an independent replica; it is not a partition of a larger job. Changing NCASC manually can change the cutoff and therefore the physics. Runs with different nuclei, excitation energies, spin inputs or normalizations must not be combined as if they were interchangeable replicas.

## Isolation and results

Each job has a directory under `OUT/jobs/ID/` containing:

- `input.pace`, copied without changes;
- `input.cs4`, HTML/RTF reports and any particle output;
- `stdout.log` and `stderr.log`;
- `runtime/` with its own executable, database, `lisepp.ini`, and `files/` directory.

The private runtime is necessary because PACE4 uses the fixed traceback filename `files/PACE4evt.dat`. The generated configuration selects `WorkInROOT=1`, so concurrent jobs cannot overwrite that file. Runtime libraries and plugins are loaded read-only from the selected installation. Large Qt libraries are not copied for every job.

`results.json` records job status, seed, original NCASC, effective cutoff, input/executable/database hashes, timings, output hashes and errors. Scientific reproducibility checks use `.cs4`, particle files and binary traceback when enabled. HTML/RTF reports include run-specific dates and paths, so their complete byte streams are not expected to match across output directories.

## Failure and cancellation

A failed or timed-out job stops the batch; active children are terminated and reaped, and jobs not yet started are marked canceled. Ctrl+C also cancels children. Logs and partial results remain available. Treat the nonzero runner exit code and per-job status as failure; the presence of some successful output files does not make the batch complete.

Keep runtime libraries and plugins unchanged during a batch. The runner is a command-line batch tool; it does not add threads to the GUI or alter the simulation algorithms.

## Validation and speed measurement

```powershell
python tests/parallel_regression.py --exe _install_MSVC/PACE4.exe
python tests/parallel_benchmark.py --exe _install_MSVC/PACE4.exe --out _build_validation/parallel_speed_new --events 40000 --rounds 3
```

The speed benchmark runs the same 12 jobs with 1, 2 and 4 workers in counterbalanced order. It verifies unchanged inputs and exact `.cs4` hashes and residue counts. End-to-end timings include Python startup, copying private runtimes, application startup, calculation and output. Report median time, range and speedup; results depend on the processor, available memory, disk, thermal state and other activity.
