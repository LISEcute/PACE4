"""Exercise a built PACE4 executable without opening windows.

Usage: python tests/cli_regression.py path/to/PACE4.exe
Requires the Qt runtime on PATH and lisecfg/AME_DB.sqlite beside the executable.
All generated inputs and outputs live in a temporary directory beside the binary.
"""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile


def input_text(mode=2, spin=10, weights=None, energy=30.125):
    lines = [f"64 0 {mode} 0 0.75 1.1 0 10.25 0 0 0",
             f"28 60 {energy} 0 {spin} 0"]
    if mode == 3:
        values = weights if weights is not None else [1] * int(spin)
        lines.append(str(len(values)) + " " + " ".join(map(str, values)))
    lines.extend(["", "", ""])
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    args = parser.parse_args()
    exe = args.executable.resolve()
    env = dict(os.environ, QT_QPA_PLATFORM="offscreen", PACE4_SEED="436")
    passed = []
    with tempfile.TemporaryDirectory(prefix="pace_cli_", dir=exe.parent) as tmp:
        root = Path(tmp)

        def run(name, text, success=True, block_rtf=False):
            path = root / (name + ".pace")
            path.write_text(text, encoding="ascii")
            if block_rtf:
                path.with_suffix(".rtf").mkdir()
            result = subprocess.run([str(exe), str(path), "-s"], cwd=root,
                                    env=env, capture_output=True, text=True,
                                    timeout=90)
            if success:
                assert result.returncode == 0, (name, result.returncode, result.stderr)
                for suffix in (".rtf", ".html", ".cs4"):
                    output = path.with_suffix(suffix)
                    assert output.is_file() and output.stat().st_size, (name, suffix)
                report = path.with_suffix(".html").read_text(encoding="utf-8")
                assert "4.38.1" in report, name + ": missing revised version"
            else:
                assert result.returncode == 1, (name, result.returncode, result.stderr)
            passed.append(name)
            return path

        first = run("fixed_spin", input_text())
        cs = first.with_suffix(".cs4").read_text()
        assert "30.125" in cs, "Loaded excitation energy changed before simulation"
        repeat = run("fixed_spin_replay", input_text())
        assert repeat.with_suffix(".cs4").read_bytes() == first.with_suffix(".cs4").read_bytes(), "Same seed/input did not replay"
        distributed = run("custom_spin", input_text(mode=3, spin=3))
        repeated = run("custom_spin_replay", input_text(mode=3, spin=3))
        assert distributed.with_suffix(".cs4").read_bytes() == repeated.with_suffix(".cs4").read_bytes(), "Distributed-spin seed replay failed"
        run("generated_spin", input_text(mode=5))
        run("spin_out_of_range", input_text(mode=5, spin=200), success=False)
        run("too_many_weights", input_text(mode=3, spin=3, weights=[1] * 201), success=False)
        run("rtf_open_failure", input_text(), success=False, block_rtf=True)
        run("truncated_input", "64 0 2\n", success=False)
    print(f"PASS: {len(passed)} CLI cases; excitation-energy preservation and seed replay verified")
    for name in passed:
        print("  " + name)


if __name__ == "__main__":
    main()
