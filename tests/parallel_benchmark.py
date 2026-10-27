"""Measure whole-job process scaling and verify scheduling-independent physics output."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import sys
import time

from physics_benchmark import CASES, input_text, residual_counts

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--events',type=int,default=40000)
    parser.add_argument('--rounds',type=int,choices=[1,2,3],default=3)
    args=parser.parse_args()
    repo=Path(__file__).resolve().parents[1]
    root=args.out.resolve()
    if root.exists(): raise SystemExit('Use a new benchmark output directory.')
    if not 1<=args.events<=1000000: raise SystemExit('Events must be between 1 and 1000000.')
    root.mkdir(parents=True)
    fixtures=root/'fixtures'; fixtures.mkdir()
    jobs=[]
    for fyrst in [0,1]:
        for case in CASES:
            path=fixtures/(case['name']+f'_fyrst{fyrst}.pace')
            path.write_text(input_text(case,args.events,fyrst),encoding='ascii')
            for replica in [0,1]:
                jobs.append({'id':case['name']+f'_f{fyrst}_r{replica}',
                             'input':path.name,'seed':436+104729*len(jobs)})
    manifest=fixtures/'jobs.json'
    manifest.write_text(json.dumps({'version':1,'jobs':jobs},indent=2),encoding='utf-8')
    record={'schema_version':1,'python':sys.version,'platform':platform.platform(),
            'logical_cpus':os.cpu_count(),'executable':str(args.exe.resolve()),
            'exe_sha256':digest(args.exe),'runner_sha256':digest(repo/'tools/pace_parallel.py'),
            'jobs_per_run':len(jobs),'events_per_job':args.events,
            'effective_probability_cutoff':1/max(500.,3*math.sqrt(args.events)+10),
            'manifest_sha256':digest(manifest),'runs':[],
            'scope':'Whole independent jobs; same input bytes and seeds; no cascade splitting.',
            'timing':'End-to-end wall time includes Python startup, runtime copying, PACE4 startup, simulation and output; no other benchmark launched concurrently.'}
    baseline=None
    orders=[[1,2,4],[4,1,2],[2,4,1]]
    for round_index,order in enumerate(orders[:args.rounds],1):
        for workers in order:
            destination=root/f'round{round_index}_workers{workers}'
            command=[sys.executable,str(repo/'tools/pace_parallel.py'),'--exe',str(args.exe.resolve()),
                     '--manifest',str(manifest),'--out',str(destination),'--jobs',str(workers),'--timeout','300']
            print(f'START round={round_index} workers={workers}',flush=True)
            started=time.perf_counter()
            completed=subprocess.run(command,capture_output=True,timeout=1800)
            elapsed=time.perf_counter()-started
            (root/f'round{round_index}_workers{workers}.log').write_bytes(completed.stdout+completed.stderr)
            if completed.returncode: raise RuntimeError(f'Runner exit {completed.returncode}; inspect {destination} and log')
            run_manifest=json.loads((destination/'results.json').read_text(encoding='utf-8'))
            assert run_manifest['status']=='succeeded'
            declared={j['id']:j for j in run_manifest['jobs']}
            signatures={}
            for job in jobs:
                folder=destination/'jobs'/job['id']
                copied=folder/'input.pace'
                assert copied.read_bytes()==(fixtures/job['input']).read_bytes(), 'Input or NCASC changed'
                report=(folder/'input.html').read_text(encoding='utf-8')
                reported_cutoff=float(re.search(r'probability discriminator of program set to ([0-9.eE+-]+)',report).group(1))
                expected_cutoff=record['effective_probability_cutoff']
                assert math.isclose(reported_cutoff,expected_cutoff,rel_tol=5e-5), 'Reported cutoff differs from unchanged NCASC'
                assert declared[job['id']]['ncasc']==args.events
                assert declared[job['id']]['seed']==job['seed']
                assert declared[job['id']]['effective_discriminator']==expected_cutoff
                counts=residual_counts(report)
                assert sum(counts.values())==args.events, (job['id'],counts)
                signatures[job['id']]={'input_sha256':digest(copied),'cs4_sha256':digest(folder/'input.cs4'),
                                       'residual_counts':counts,'seed':job['seed'],'reported_cutoff':reported_cutoff}
            if baseline is None: baseline=signatures
            assert signatures==baseline, f'Physics outputs changed with scheduling in {destination}'
            record['runs'].append({'round':round_index,'workers':workers,'wall_seconds':elapsed,
                                   'job_outputs_identical':True,'folder':str(destination)})
            (root/'timing.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
            print(f'PASS round={round_index} workers={workers}: {elapsed:.3f}s; all {len(jobs)} inputs and CS4 hashes identical',flush=True)
    serial=statistics.median(r['wall_seconds'] for r in record['runs'] if r['workers']==1)
    record['summary']=[]
    for workers in [1,2,4]:
        samples=[r['wall_seconds'] for r in record['runs'] if r['workers']==workers]
        med=statistics.median(samples)
        record['summary'].append({'workers':workers,'samples_seconds':samples,'median_seconds':med,
                                  'min_seconds':min(samples),'max_seconds':max(samples),
                                  'speedup':serial/med,'parallel_efficiency':serial/med/workers})
    (root/'timing.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    (root/'physics_signatures.json').write_text(json.dumps(baseline,indent=2),encoding='utf-8')
    print(json.dumps(record['summary'],indent=2),flush=True)

if __name__=='__main__': main()
