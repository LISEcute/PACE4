"""Test the independent-process runner, optionally including a deployed PACE4 binary.

python tests/parallel_regression.py
python tests/parallel_regression.py --exe _install_MSVC/PACE4.exe
The real-app check runs trace-enabled IDIST=1/2 jobs both serially and concurrently.
"""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('pace_parallel', ROOT / 'tools/pace_parallel.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
REAL_POPEN = subprocess.Popen
REAL_EXE = None

FAKE_WORKER = r'''
import os,pathlib,sys,time
p=pathlib.Path(sys.argv[1]);data=p.read_bytes();control=data.decode().splitlines()[1]
if control=='sleep':time.sleep(30)
if control=='fail':sys.exit(7)
if control=='missing':sys.exit(0)
assert 'PACE_BENCHMARK_SPINS' not in os.environ
payload=data+b'\nseed='+os.environ['PACE4_SEED'].encode()
for suffix in ('.cs4','.html','.rtf'):p.with_suffix(suffix).write_bytes(payload)
pathlib.Path('runtime/files/PACE4evt.dat').write_bytes(payload)
'''


def scientific_hashes(result):
    return {j['id']: {name: metadata['sha256'] for name, metadata in j['outputs'].items()
                     if name.endswith(('.cs4', '.particles', '.dat'))} for j in result['jobs']}


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='pace_parallel_test_')
        def cleanup():
            # Windows can briefly retain file locks after a terminated child exits.
            for attempt in range(30):
                try:
                    self.temporary.cleanup()
                    return
                except PermissionError:
                    if attempt == 29: raise
                    time.sleep(.1)
        self.addCleanup(cleanup)
        self.root = Path(self.temporary.name)
        self.database = self.root / 'lisecfg'
        self.database.mkdir()
        (self.database / 'AME_DB.sqlite').write_bytes(b'independent fake mass data')
        self.worker = self.root / 'worker.py'
        self.worker.write_text(FAKE_WORKER)
        self.processes = []
        self.sequence = 0

    def manifest(self, controls=('ok', 'ok', 'ok'), events=100000):
        jobs = []
        for i, control in enumerate(controls):
            source = self.root / f'source-{i}.pace'
            source.write_text(f'{events} 1 2 0 -1 1 0 10 0 1 0\n{control}\n')
            jobs.append(dict(id=f'job-{i}', input=source.name, seed=436+i*104729))
        path = self.root / 'manifest.json'
        path.write_text(json.dumps(dict(version=1, jobs=jobs)))
        return path

    def spawn(self, command, **options):
        # Real, isolated child processes; only substitute the native PACE payload.
        p = REAL_POPEN([sys.executable, str(self.worker), *command[1:]], **options)
        self.processes.append(p)
        return p

    def run_fake(self, manifest, workers=1, timeout=None, spawn=None):
        self.sequence += 1
        out = self.root / f'run-{self.sequence}'
        with patch.object(runner.subprocess, 'Popen', side_effect=spawn or self.spawn):
            result = runner.run_suite(Path(sys.executable), manifest, out, workers, timeout, self.database)
        self.assertTrue(all(p.poll() is not None for p in self.processes), 'Every child must be reaped')
        return result

    def test_serial_parallel_and_manifest_order_reproduce(self):
        manifest = self.manifest()
        with patch.dict(os.environ, {'PACE4_SEED': '999', 'PACE_BENCHMARK_SPINS': 'must-not-be-used'}):
            c1, serial = self.run_fake(manifest, 1)
            c2, parallel = self.run_fake(manifest, 3)
        document = json.loads(manifest.read_text())
        document['jobs'].reverse()
        manifest.write_text(json.dumps(document))
        c3, reordered = self.run_fake(manifest, 2)
        self.assertEqual((c1,c2,c3), (0,0,0))
        self.assertEqual(scientific_hashes(serial), scientific_hashes(parallel))
        self.assertEqual(scientific_hashes(serial), scientific_hashes(reordered))
        for job in parallel['jobs']:
            folder=Path(job['folder'])
            self.assertEqual((folder/'input.pace').read_bytes(), Path(job['input_source']).read_bytes())
            self.assertEqual(job['ncasc'], 100000)
            self.assertEqual(job['effective_discriminator'], 1/(3*math.sqrt(100000)+10))
            self.assertIn('WorkInROOT=1', (folder/'runtime/lisepp.ini').read_text())
            self.assertEqual((folder/'runtime/lisecfg/AME_DB.sqlite').read_bytes(), b'independent fake mass data')
        self.assertEqual(len({j['folder'] for j in parallel['jobs']}),3)

    def test_reject_manifest_errors_without_creating_outputs(self):
        manifest=self.manifest()
        original=json.loads(manifest.read_text())
        bad=[]
        for key,value in [('id','../escape'),('id','CON'),('seed',0),('seed',True),('seed',2147483647)]:
            d=json.loads(json.dumps(original));d['jobs'][0][key]=value;bad.append(d)
        d=json.loads(json.dumps(original));d['jobs'][1]['id']='JOB-0';bad.append(d)
        d=json.loads(json.dumps(original));d['jobs'][1]['seed']=d['jobs'][0]['seed'];bad.append(d)
        d=json.loads(json.dumps(original));d['jobs'][0]['seed']=runner.RNG_MASK;d['jobs'][1]['seed']=runner.RNG_MASK^1;bad.append(d)
        d=json.loads(json.dumps(original));d['jobs'][0]['ncasc']=5;bad.append(d)
        d=json.loads(json.dumps(original));d['split']=2;bad.append(d)
        for document in bad:
            manifest.write_text(json.dumps(document))
            with self.assertRaises(ValueError):
                runner.run_suite(Path(sys.executable),manifest,self.root/'forbidden',data_dir=self.database)
            self.assertFalse((self.root/'forbidden').exists())

    def test_existing_output_refused_and_missing_database(self):
        manifest=self.manifest()
        existing=self.root/'existing';existing.mkdir();(existing/'keep.txt').write_text('keep')
        with self.assertRaises(FileExistsError):
            runner.run_suite(Path(sys.executable),manifest,existing,data_dir=self.database)
        self.assertEqual((existing/'keep.txt').read_text(),'keep')
        with self.assertRaises(ValueError):
            runner.run_suite(Path(sys.executable),manifest,self.root/'new',data_dir=self.root/'absent')
        self.assertFalse((self.root/'new').exists())

    def test_failure_cancels_active_and_pending_children(self):
        code,result=self.run_fake(self.manifest(('fail','sleep','ok')),2)
        self.assertEqual(code,1)
        self.assertEqual([j['status'] for j in result['jobs']],['failed','canceled','canceled'])
        self.assertEqual(result['jobs'][0]['returncode'],7)

    def test_timeout_terminates_children(self):
        code,result=self.run_fake(self.manifest(('sleep','sleep','ok')),2,.15)
        self.assertEqual(code,1)
        self.assertIn('timed_out',[j['status'] for j in result['jobs']])
        self.assertEqual(result['jobs'][-1]['status'],'canceled')

    def test_success_without_outputs_is_failure(self):
        code,result=self.run_fake(self.manifest(('missing',)))
        self.assertEqual(code,1)
        self.assertEqual(result['jobs'][0]['status'],'failed')
        self.assertIn('missing/empty',result['jobs'][0]['error'])

    def test_launch_error_stops_already_running_child(self):
        def spawn(command,**options):
            if self.processes: raise OSError('simulated launch failure')
            return self.spawn(command,**options)
        code,result=self.run_fake(self.manifest(('sleep','ok','ok')),2,spawn=spawn)
        self.assertEqual(code,1)
        self.assertEqual([j['status'] for j in result['jobs']],['canceled','failed','canceled'])

    def test_interrupt_stops_children_and_records_summary(self):
        def interrupt():
            deadline=time.monotonic()+5
            while len(self.processes)<2 and time.monotonic()<deadline:time.sleep(.01)
            signal.raise_signal(signal.SIGINT)
        thread=threading.Thread(target=interrupt)
        thread.start()
        try:
            code,result=self.run_fake(self.manifest(('sleep','sleep','ok')),2)
        finally:
            thread.join(timeout=6)
        self.assertEqual(code,130)
        self.assertEqual(result['status'],'canceled')
        self.assertTrue(all(j['status']=='canceled' for j in result['jobs']))
        persisted=json.loads((Path(result['jobs'][0]['folder']).parents[1]/'results.json').read_text())
        self.assertEqual(persisted['status'],'canceled')

    def test_header_and_worker_validation(self):
        manifest=self.manifest(events=0)
        with self.assertRaises(ValueError):runner.read_manifest(manifest)
        manifest=self.manifest(events=1000001)
        with self.assertRaises(ValueError):runner.read_manifest(manifest)
        manifest=self.manifest()
        for kwargs in ({'workers':0},{'workers':True},{'timeout':0},{'timeout':float('inf')}):
            with self.assertRaises(ValueError):
                runner.run_suite(Path(sys.executable),manifest,self.root/'invalid',data_dir=self.database,**kwargs)

    def test_real_pace_trace_enabled_serial_parallel(self):
        if REAL_EXE is None: self.skipTest('Supply --exe to test the deployed PACE4 binary')
        jobs=[]
        for idist in (0,1,2):
            source=self.root/f'trace-{idist}.pace'
            source.write_text(f'128 {idist} 5 1 -1 1 0 10 0 1 31\n28 60 30.125 0 20 0\n\n\n\n')
            jobs.append(dict(id=f'trace-{idist}',input=source.name,seed=104729+idist*130363))
        manifest=self.root/'real.json';manifest.write_text(json.dumps(dict(version=1,jobs=jobs)))
        c1,serial=runner.run_suite(REAL_EXE,manifest,self.root/'real-serial',1,90)
        c2,parallel=runner.run_suite(REAL_EXE,manifest,self.root/'real-parallel',3,90)
        self.assertEqual((c1,c2),(0,0), (serial,parallel))
        self.assertEqual(scientific_hashes(serial),scientific_hashes(parallel))
        traces=[]
        for job in parallel['jobs'][1:]:
            event=Path(job['folder'])/'runtime/files/PACE4evt.dat'
            self.assertTrue(event.is_file() and event.stat().st_size>0)
            traces.append(event)
        self.assertNotEqual(traces[0],traces[1])
        malformed=self.root/'malformed.pace'
        malformed.write_text('128 0 2 0 -1 1 0 10 0 1 0\n')
        bad_manifest=self.root/'real-failure.json'
        bad_manifest.write_text(json.dumps(dict(version=1,jobs=[
            dict(id='a-malformed',input=malformed.name,seed=7),
            dict(id='b-pending',input=jobs[0]['input'],seed=8)])))
        failed,summary=runner.run_suite(REAL_EXE,bad_manifest,self.root/'real-failure',1,90)
        self.assertEqual(failed,1)
        self.assertEqual([j['status'] for j in summary['jobs']],['failed','canceled'])
        print('Real PACE4: IDIST0/1/2 scientific and binary-event hashes match serial/parallel; trace paths isolated; malformed input fails safely.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe',type=Path)
    args,remaining=parser.parse_known_args()
    if args.exe:
        REAL_EXE=args.exe.resolve()
    unittest.main(argv=[sys.argv[0],*remaining])
