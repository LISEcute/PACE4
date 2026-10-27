"""Build/run isolated GNU-Fortran PACE3 references; never modify historical sources.

Portability adaptations: VAX OPTIONS removal; -fno-automatic static storage;
replacement of VAX RAN(seed) by the same portable masked Park-Miller recurrence
as PACE4 (GNU RAN does NOT update seed and otherwise returns a constant);
seed environment input and read-only event instrumentation. The checked build
also avoids a zero-weight out-of-bounds FISROT interpolation neighbor. Compile
an unchecked companion without this last repair to measure its numerical effect.
"""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

SEEDS = [436, 104729, 130363, 155921, 196613]
CASES = [dict(name="ni60_fixed_j10", mode=2, z=28, a=60, energy=30.125, spin=10),
         dict(name="ni60_equal_j0_10_20", mode=3, z=28, a=60, energy=30.125, spin=20),
         dict(name="ni60_generated_j0_20", mode=5, z=28, a=60, energy=30.125, spin=20)]


def reference_source(original, repair):
    source = re.sub(r'^      options /g_floating.*$',
                    'C     VAX OPTIONS /G_FLOATING omitted for GNU Fortran',
                    original, flags=re.M|re.I)
    old = """      function ranf()
      integer*4 xran
      data xran/30031961/
      ranf=ran(xran)
      return
      end"""
    new = """      function ranf()
C     Portable substitute for VAX RAN(seed), REAL*4 result retained.
      integer*8 state,k,mask,im
      integer ios
      logical first
      character*32 seedtext
      parameter(mask=123459876_8,im=2147483647_8)
      data state/30031961_8/,first/.true./
      if(first) then
      call get_environment_variable('PACE4_SEED',seedtext)
      if(len_trim(seedtext).gt.0) read(seedtext,*,iostat=ios)state
      if(state.le.0.or.state.ge.im) stop 2
      write(9,*) 'Portable masked Park-Miller seed: ',state
      if(ieor(state,mask).eq.0.or.ieor(state,mask).eq.im)
     1 state=ieor(1_8,mask)
      first=.false.
      endif
      state=ieor(state,mask)
      k=state/127773_8
      state=16807_8*(state-k*127773_8)-2836_8*k
      if(state.lt.0)state=state+im
      ranf=real(dble(state)/dble(im))
      state=ieor(state,mask)
      return
      end"""
    assert old in source
    source = source.replace(old, new)
    initial = """      OPEN(71,FILE='initial_spins.dat',STATUS='REPLACE')
      DO 9871 K=1,NCASC
 9871 WRITE(71,*) K,JCS(K),MJCS(K),ERGCS(K)
      CLOSE(71)
"""
    needle = '   17     CONTINUE\n      GO TO  100'
    assert needle in source
    source = source.replace(needle, '   17     CONTINUE\n'+initial+'      GO TO  100')
    source = source.replace(' 1308 JCS(L)=MAXC', ' 1308 JCS(L)=MAXC\n'+initial.replace('9871','9872'))
    source = source.replace('      CALL PRODCT(II,IDIST,IZMAX,INMAX,NCASC)', """      OPEN(72,FILE='final_events.dat',STATUS='REPLACE')
      DO 9873 K=1,NCASC
 9873 WRITE(72,*) K,IZCS(K),INCS(K),JCS(K),ERGCS(K)
      CLOSE(72)
      OPEN(73,FILE='normalization.dat',STATUS='REPLACE')
      WRITE(73,*) SIGMA,NCASC
      CLOSE(73)
      CALL PRODCT(II,IDIST,IZMAX,INMAX,NCASC)""")
    source = source.replace("status='new'", "status='replace'")
    if repair:
        old = '      B2=(X2B(IX+1,IY+1 )-X2B(IX,IY+1))*DX+X2B(IX,IY+1)'
        assert old in source
        source = source.replace(old, '      B2=B1\n      IF(IY.LT.11) THEN\n'+old+'\n      ENDIF')
    return source


def input_text(case, events, fyrst=0):
    lines = [f"{events} 0 {case['mode']} {fyrst} -1 1 0 10 0 1",
             f"{case['z']:5d}{case['a']:5d}{case['energy']:10.5f}{0.:10.5f}{case['spin']:10.5f}{0:5d}"]
    if case['mode'] == 3:
        weights = [100/3 if j in (0,10,20) else 0. for j in range(21)]
        lines.append(f'{len(weights):5d}'+' '*5+''.join(f'{v:10.5f}' for v in weights[:7]))
        for start in range(7,len(weights),8):
            lines.append(''.join(f'{v:10.5f}' for v in weights[start:start+8]))
    # Three default optical potentials, default EM strengths, no discrete levels.
    return '\n'.join(lines+['','','','','    0    0    0'])+'\n'


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--gfortran',type=Path,default=Path('C:/Qt/Tools/mingw1120_64/bin/gfortran.exe'))
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--events',type=int,default=5000)
    ap.add_argument('--replicas',type=int,default=5,choices=range(1,6))
    ap.add_argument('--unchecked',action='store_true')
    ap.add_argument('--seed-offset',type=int,default=1000003)
    ap.add_argument('--fyrst',type=float,default=0.)
    args=ap.parse_args()
    assert 0<args.events<=5016
    repo=Path(__file__).resolve().parents[1]
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    original_path=repo/'g_PACE/fortran/PACE3.for'
    original=original_path.read_text()
    source=out/'pace3_benchmark.f';source.write_text(reference_source(original,not args.unchecked))
    exe=out/'pace3_benchmark.exe'
    env=dict(os.environ);env['PATH']=str(args.gfortran.parent)+';'+env['PATH']
    flags=['-std=legacy','-fdec','-fallow-argument-mismatch','-ffixed-line-length-72','-fno-automatic','-O0']
    if not args.unchecked: flags+=['-fcheck=bounds']
    command=[str(args.gfortran)]+flags+['-o',str(exe),str(source)]
    build=subprocess.run(command,env=env,capture_output=True)
    (out/'compile.log').write_bytes(build.stdout+build.stderr)
    if build.returncode: raise RuntimeError('Reference build failed; see compile.log')
    suite=dict(executable=str(exe),executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
               original_sha256=hashlib.sha256(original_path.read_bytes()).hexdigest(),
               adapted_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
               compile_command=command,reference='PACE3.for with documented GNU portability adaptations',
               bounds_repair=not args.unchecked,seed_offset=args.seed_offset,events_per_run=args.events,fyrst=args.fyrst,
               noshl=1,facla=10,aratio=1,barfac=-1,runs=[])
    for case in CASES:
        for seed in SEEDS[:args.replicas]:
            folder=out/case['name']/str(seed);folder.mkdir(parents=True,exist_ok=True)
            (folder/'pace3in.dat').write_text(input_text(case,args.events,args.fyrst))
            env['PACE4_SEED']=str(seed+args.seed_offset)
            start=time.monotonic()
            result=subprocess.run([str(exe)],cwd=folder,env=env,capture_output=True,timeout=180)
            duration=time.monotonic()-start
            (folder/'stdout.log').write_bytes(result.stdout);(folder/'stderr.log').write_bytes(result.stderr)
            if result.returncode or not (folder/'final_events.dat').exists():
                raise RuntimeError(f'{case["name"]}/{seed}: failed; inspect {folder}')
            spins=Counter(int(line.split()[1])-1 for line in (folder/'initial_spins.dat').read_text().splitlines())
            residuals=Counter()
            for line in (folder/'final_events.dat').read_text().splitlines():
                fields=line.split();z,n=int(fields[1]),int(fields[2])
                residuals['fission' if z<0 else f'{z},{n}']+=1
            sigma=float((folder/'normalization.dat').read_text().split()[0])
            assert sum(spins.values())==args.events and sum(residuals.values())==args.events
            suite['runs'].append(dict(case=case,seed=seed,actual_seed=seed+args.seed_offset,events=args.events,seconds=duration,
                                      sigma_input_mb=sigma,initial_spin_counts=dict(spins),
                                      residual_counts=dict(residuals),folder=str(folder)))
            (out/'results.json').write_text(json.dumps(suite,indent=2))
            print(f'{case["name"]} seed={seed}: {duration:.2f}s',flush=True)

if __name__=='__main__': main()
