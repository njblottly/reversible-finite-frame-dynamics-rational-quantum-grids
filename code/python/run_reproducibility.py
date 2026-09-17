#!/usr/bin/env python3
"""Save a dated local record for both RaQM Python verification scripts.

Place this file alongside verify_gate_composition.py and verify_nominal_gates.py.
Run: python3 run_reproducibility.py
No packages, network access, administrator privileges, or shell piping needed.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime,timezone
from pathlib import Path

def now():return datetime.now(timezone.utc).isoformat()
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--script-dir',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--records-dir',type=Path,default=None)
    args=parser.parse_args()
    base=args.script_dir.resolve()
    names=('verify_gate_composition.py','verify_nominal_gates.py')
    missing=[n for n in names if not (base/n).is_file()]
    if missing:
        parser.error('Put these files in --script-dir first: '+', '.join(missing))
    root=args.records_dir.resolve() if args.records_dir else base/'reproducibility'
    out=root/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out.mkdir(parents=True,exist_ok=False)
    record={'started_utc':now(),'python_version':sys.version,
        'python_executable':sys.executable,'platform':platform.platform(),
        'machine':platform.machine(),'macos_version':platform.mac_ver()[0],
        'runner_sha256':sha(Path(__file__).resolve()),'script_directory':str(base),
        'runs':[],'mathematica':'Not executed by this Python runner.'}
    for name in names:
        script=base/name
        output=out/(script.stem+'_results.json')
        command=[sys.executable,str(script),'--output',str(output)]
        entry={'script':name,'sha256':sha(script),'started_utc':now(),'command':command}
        try:
            completed=subprocess.run(command,cwd=base,capture_output=True,text=True,check=False)
            entry.update({'finished_utc':now(),'returncode':completed.returncode})
            (out/(script.stem+'_stdout.txt')).write_text(completed.stdout,encoding='utf-8')
            (out/(script.stem+'_stderr.txt')).write_text(completed.stderr,encoding='utf-8')
            entry['status']='passed' if completed.returncode==0 else 'failed'
            if output.is_file():entry['result_sha256']=sha(output)
            print(f"{name}: {entry['status']}")
            if completed.returncode:print(completed.stderr)
        except OSError as exc:
            entry.update({'finished_utc':now(),'status':'failed','exception':str(exc)})
            print(f'{name}: failed: {exc}')
        record['runs'].append(entry)
        (out/'run_record.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    evidence=base/'mathematica_verification_record.json'
    if evidence.is_file():
        (out/evidence.name).write_bytes(evidence.read_bytes())
        record['prior_mathematica_evidence_sha256']=sha(evidence)
    record['finished_utc']=now()
    record['status']='passed' if all(r['status']=='passed' for r in record['runs']) else 'failed'
    (out/'run_record.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print('Local record:',out)
    print('Overall status:',record['status'])
    raise SystemExit(0 if record['status']=='passed' else 1)

if __name__=='__main__':main()
