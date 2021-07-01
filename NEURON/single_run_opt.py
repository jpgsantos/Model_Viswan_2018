# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 17:48:41 2016

@author: daniel
"""

from neuron import h
import d1msn as msn
import EGF_experiment as pe
import pickle
import parameters as p
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
#import single_sec_cell as se

# --- 1. Create a cell and other useful stuff

with open('D1_71bestFit_updRheob.pkl', 'rb') as f:
    model_sets  = pickle.load(f, encoding="latin1")

cell_ID = 70
variables = model_sets[cell_ID]['variables'] 
cell = msn.MSN(variables = variables) 

#cell = se.SingleSecCell()

cell.somalist[0].insert('Viswan_2018_opt')
cell.somalist[0].EGF_start_Viswan_2018_opt = p.EGF_input_start
cell.somalist[0].EGF_length_Viswan_2018_opt = p.EGF_input_length
cell.somalist[0].EGF_level_Viswan_2018_opt = p.EGF_level
cell.somalist[0].EGF_steepness_Viswan_2018_opt = p.EGF_steepness
# Set the pointer to MAPK_P
#for sec in h.allsec():
#    if 'axon' in sec.name():    
#        for seg in sec:
#            seg.total_MAPK_Im = p.total_MAPK
#            h.setpointer(cell.somalist[0](0.5)._ref_MAPK_p_Viswan_2018, 'MAPK_P', seg.Im)

#import pandas as pd
#df = pd.read_csv('E0_pERK.csv')
#MAPK_P = df['MAPK_P'].values.tolist()
#MAPK_P_P = df['MAPK_P_P'].values.tolist()
#pERK = np.asarray(MAPK_P) + np.asarray(MAPK_P_P)
#plt.plot(pERK)

#h('MAPK_P = 0.00456402')
#h('MAPK_P = 0.0')
#h('MAPK_P_P = 0.00456402')
for sec in h.allsec():
    if not 'axon' in sec.name():    
        for seg in sec:
            seg.total_MAPK_kaf = p.total_MAPK
            h.setpointer(cell.somalist[0](0.5)._ref_MAPK_p_Viswan_2018_opt, 'MAPK_P', seg.kaf)
            h.setpointer(cell.somalist[0](0.5)._ref_MAPK_p_p_Viswan_2018_opt, 'MAPK_P_P', seg.kaf)
#            h.setpointer(h._ref_MAPK_P, 'MAPK_P', seg.kaf)
#            h.setpointer(h._ref_MAPK_P_P, 'MAPK_P_P', seg.kaf)
#            h.setpointer(MAPK_P, 'MAPK_P', seg.kaf)
#            h.setpointer(MAPK_P_P, 'MAPK_P_P', seg.kaf)

#
gkaf = h.Vector()
gkaf.record(cell.dendlist[10](0.5)._ref_gbarmod_kaf, p.record_step)
# --- 2. Insert stimulation to cell

#independent_dends = [3, 5, 8, 12, 15, 22, 26, 35, 41, 47, 53, 57]
dend_record_list = [3]
plateau_cluster_list = [3]

ex = pe.EGF_Experiment('record_ca_opt', cell)           
#ex.insert_synapses('MSN')
ex.set_up_recording()  
ex.simulate()
ex.plot_results()
plt.plot(ex.tout, gkaf)
plt.show()