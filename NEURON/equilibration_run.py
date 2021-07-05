# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 17:48:41 2016

@author: daniel
"""

import single_sec_cell as se
import EGF_experiment as te
import json

opt_flag = True
simtime = 5000e3
cell = se.SingleSecCell(opt_flag = opt_flag)
if not opt_flag:
    cell.somalist[0].EGF_level_Findsim = 0.0
    cell.somalist[0].EGF_steepness_Findsim = 1e-5
else:
    cell.somalist[0].EGF_level_Findsim_opt = 0.0
    cell.somalist[0].EGF_steepness_Findsim_opt = 1e-5
ex = te.EGF_Experiment('single_sec', cell)           
ex.set_up_recording()
ex.simulate(simtime = simtime)
ex.plot_results()

steady_states = []
for sp in ex.species:
    print(sp.to_python()[-1])
    steady_states.append(sp.to_python()[-1])

result = {'steady_states': steady_states}
to_save = json.dumps(result)
if not opt_flag:
    filename = 'steady_states.dat'
else:
    filename = 'steady_states_opt.dat'
with open(filename, 'w', encoding = 'utf-8') as f:
    json.dump(to_save, f)
    

