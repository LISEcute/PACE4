"""Replay default high-spin cases and earlier low-spin cases across two builds."""
import argparse, hashlib, json, sys
from pathlib import Path
repo=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(repo/'tools'))
from pace_parallel import run_suite
from physics_benchmark import CASES, input_text, residual_counts

def main():
 p=argparse.ArgumentParser()
 p.add_argument('--exe',type=Path,required=True)
 p.add_argument('--out',type=Path,required=True)
 p.add_argument('--compare',type=Path)
 a=p.parse_args(); out=a.out.resolve(); out.mkdir(parents=True,exist_ok=False)
 optical='47 -0.267 -0.002 1.264 0.660 0 9.520 -0.053 0 1.332 0.480 0 250 0 1\n'
 jobs=[]
 for fyrst in [0,1]:
  for seed in [436,104729]:
   name=f'default_fyrst{fyrst}_seed{seed}'
   data=f'1000 1 1 {fyrst} 0 1 0 10 0 0 0\n20 48 50 124 0 0 0\n300 0 0 2 0 0 0 3\n'+optical*4+'\n\n'
   # Unique external seeds across rows; labels record the actual chosen value.
   actual=seed+fyrst*1000003
   (out/f'{name}.pace').write_text(data,encoding='ascii')
   jobs.append(dict(id=name,input=f'{name}.pace',seed=actual))
 for i,case in enumerate(CASES):
  name=case['name']; (out/f'{name}.pace').write_text(input_text(case,5000),encoding='ascii')
  jobs.append(dict(id=name,input=f'{name}.pace',seed=196613+i*104729))
 manifest=out/'jobs.json'; manifest.write_text(json.dumps(dict(version=1,jobs=jobs)),encoding='utf-8')
 code,result=run_suite(a.exe,manifest,out/'runs',2,180)
 assert code==0, result.get('error','Batch failed')
 evidence={}
 for j in result['jobs']:
  folder=Path(j['folder']); report=(folder/'input.html').read_text(encoding='utf-8')
  assert '"FIND" limits error' not in report,j['id']+': FIND warning'
  assert 'LEVEL DENSITY EXCEEDED' not in report,j['id']+': table overflow'
  counts=residual_counts(report)
  assert sum(counts.values())==j['ncasc'],(j['id'],counts)
  evidence[j['id']]=dict(seed=j['seed'],counts=counts,cs4_sha256=hashlib.sha256((folder/'input.cs4').read_bytes()).hexdigest())
 if a.compare:
  other=json.loads(a.compare.read_text())
  assert evidence==other,'Cross-build scientific results differ'
 (out/'scientific_results.json').write_text(json.dumps(evidence,indent=2),encoding='utf-8')
 print(f'PASS {len(jobs)} cases: no FIND/capacity warnings, event totals conserved'+('; identical cross-build results' if a.compare else ''))

if __name__=='__main__': main()
