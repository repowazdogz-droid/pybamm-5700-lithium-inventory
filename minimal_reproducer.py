"""Wrapper reproducing the 20/20, rtol=atol=1e-8 configurations quoted in the comment.

Run in an environment containing the desired PyBaMM revision and its declared
solver dependencies. Requires measurement.py and referee.py beside this file.
No plotting or GitHub access. Writes only to the specified local directory.
"""
import argparse
import json
import pathlib
import pybamm
from measurement import prepare, extract
from referee import controls

parser=argparse.ArgumentParser()
parser.add_argument('--output',type=pathlib.Path,default=pathlib.Path('reviewer-output'))
args=parser.parse_args()
args.output.mkdir(parents=True,exist_ok=True)
(args.output/'controls.json').write_text(json.dumps(controls(),indent=2))
for case in ['basic_variable','basic_constant','half']:
    simulation,p,domains,electrodes,_=prepare(case,20,20,1e-8)
    solution=simulation.solve(initial_soc=0.45)
    result=extract(simulation,solution,p,domains,electrodes,case,args.output/case)
    print(json.dumps({'case':case,'source':pybamm.__file__,'max_abs_residual_mol':result['max_abs_residual_mol'],'relative_max':result['relative_max'],'final_residual_mol':result['final_residual_mol']}))
