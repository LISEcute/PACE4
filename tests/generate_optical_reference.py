"""Reproduce optical golden values from unchanged historical Fortran routines.

Only the standalone driver is new. Promoting default REAL to 8 bytes lets the
reference and C++ use the same precision. No changes to TLOM or OPTPOT are made.
Example: python tests/generate_optical_reference.py --out _build_validation/optical_reference
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

DRIVER = """      PROGRAM OPTICAL_REFERENCE
      DIMENSION TL(150)
      CALL OPTPOT(1,28,32)
      CALL TLOM(0,1,1.,0.,59.,27.,10.,LL,TL)
      WRITE(*,*) 'neutron',LL,(TL(I),I=1,6)
      CALL OPTPOT(2,28,32)
      CALL TLOM(0,2,1.,1.,59.,27.,10.,LL,TL)
      WRITE(*,*) 'proton',LL,(TL(I),I=1,6)
      CALL OPTPOT(3,28,32)
      CALL TLOM(0,3,4.,2.,56.,26.,15.,LL,TL)
      WRITE(*,*) 'alpha',LL,(TL(I),I=1,6)
      END
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gfortran', type=Path,
                        default=Path('C:/Qt/Tools/mingw1120_64/bin/gfortran.exe'))
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    original = Path(__file__).resolve().parents[1] / 'g_PACE/fortran/PACE3.for'
    source = original.read_text()
    optical = source[source.index('      SUBROUTINE TLOM('):source.index('      SUBROUTINE MOMENT(')]
    defaults = source[source.index('      SUBROUTINE OPTPOT('):source.index('      SUBROUTINE FACLOG')]
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    adapted = out / 'optical_reference.f'
    adapted.write_text(DRIVER + optical + defaults)
    executable = out / 'optical_reference.exe'
    command = [str(args.gfortran), '-std=legacy', '-ffixed-line-length-72',
               '-fno-automatic', '-fdefault-real-8', '-fdefault-double-8',
               '-fcheck=bounds', '-O0', '-o', str(executable), str(adapted)]
    env = dict(os.environ)
    env['PATH'] = str(args.gfortran.parent) + os.pathsep + env.get('PATH', '')
    build = subprocess.run(command, env=env, capture_output=True)
    (out / 'compile.log').write_bytes(build.stdout + build.stderr)
    build.check_returncode()
    run = subprocess.run([str(executable)], env=env, capture_output=True, timeout=60)
    (out / 'coefficients.txt').write_bytes(run.stdout + run.stderr)
    run.check_returncode()
    provenance = dict(original=str(original), original_sha256=hashlib.sha256(original.read_bytes()).hexdigest(),
                      extraction='Unmodified TLOM and OPTPOT; new standalone driver; default REAL promoted to 8 bytes',
                      generated_source_sha256=hashlib.sha256(adapted.read_bytes()).hexdigest(),
                      executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
                      compile_command=command)
    (out / 'provenance.json').write_text(json.dumps(provenance, indent=2))
    print(run.stdout.decode())


if __name__ == '__main__':
    main()
