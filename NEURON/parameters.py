# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 20:39:04 2016

@author: daniel
"""
#from math import sqrt
results_directory = './results'

#-----------------------------------------------------------#
#      1. General recording and simulation parameters       #
#-----------------------------------------------------------#
step = 20.0
record_step = 5000
v_record_step = 0.1

nrn_dots_per_1ms = 1/record_step
time_to_avg_over = 20 # in seconds

simtime = 4500e3 
num_trials = 20
NUMBER_OF_PROCESSES = 7

#-----------------------------------#
#      2. Subcellular cascade       #
#             parameters            #
#-----------------------------------#

EGF_input_start = 500e3
EGF_input_length = 4000e3
EGF_level = 0.1
EGF_steepness = 0.3e-3
total_MAPK = 0.360000000000505

#-----------------------------------#
#      2. Synaptic parameters       #
#-----------------------------------#

esyn_tau = 6
isyn_tau = 6
isyn_tau_plateau = 87
e_esyn = 0
e_gaba = -60
erate = 1.0
irate = 1.0
pos = 0.5

# NMDA parameters
Mg = 1.0
alpha = 0.062
eta = 1.0/3.57

g_ramp_max = 0.000255
nmda_ampa_ratio = 1
gAMPAmax = 0.1e-3
gNMDAmax = gAMPAmax*nmda_ampa_ratio
gGABAmax = 1.5e-3
g_expsyn_max =  0.1e-3
g_inhexpsyn_max = gGABAmax

gAMPAmax_plateau = 0.5e-3 
gNMDAmax_plateau = 1.5e-3 
gGABAmax_plateau = 0.0015*1
nmda_ampa_ratio = gNMDAmax_plateau/gAMPAmax_plateau
ratio_glutamate_syn = 1.0

gAMPAmax_pf = 0.05e-3
gNMDAmax_pf = 1.5e-3 

glu_thresh1 = 0.5
glu_thresh2 = 0.5

tNMDAon = 2.76
tNMDAoff = 115.5

taur_cadyn_nmda = 160
scale_cadyn_nmda = 0.083

iclamp_amp = 1.85
iclamp_delay = 5
iclamp_periodic_delay = 100
iclamp_dur = 15
num_iclamps = 10

tau1_exp2syn = 0.1
tau2_exp2syn = 2.0

v_thresh = -50
glu_thresh = 0.06
#-----------------------------------------#
#      3. Synaptic input parameters       #
#-----------------------------------------#

plateau_syn_rate = 10
plateau_burst_start = 100
plateau_burst_end = 130
plateau_cluster_size = 10
cluster_start_pos = 0.5
cluster_end_pos = 0.5


inhibitory_syn_rate = 85.0
inhibitory_burst_start = 100
inhibitory_burst_end = 160
inhibitory_cluster_size = 1

distributed_input_rate = 1000.0/40
distributed_input_start = 300
distributed_input_end = 500
distributed_input_size = 30

deterministic_interval = 100
num_deterministic_spikes = 10

e_interval = 1.0/erate*(10**3)
i_interval = 1.0/irate*(10**3)
plateau_syn_interval = 1.0/plateau_syn_rate*(10**3)
distributed_input_interval = 1.0/distributed_input_rate*(10**3)
inhibitory_syn_interval = 1.0/inhibitory_syn_rate*(10**3)

#--------------------------------#
#      4. Spine parameters       #
#--------------------------------#
head_L = 0.5
head_diam = 0.5
neck_L = 0.5
neck_diam = 0.12
neck_Ra = 1130.0
head_Ra = 150

#-----------------------------------------------------------#
#      5. XOR problem and adaptive synapse parameters       #
#-----------------------------------------------------------#

event_times = [200, 500, 900]


#----------------------------------#
#      6. Plotting parameters      #
#----------------------------------#

#-----------------------------------#
#      6.1. For XOR experiment      #
#-----------------------------------#

#-------------------------------------------------------#
#      7. Miscellaneous and parameter dictionaries      #
#-------------------------------------------------------#

dends_per_plot = 3


params = {
        'erate' : erate,
        'irate' : irate,
        'e_esyn' : e_esyn,
        'g_expsyn_max' : g_expsyn_max,
        'g_inhexpsyn_max' : g_inhexpsyn_max,
        'erate' : erate,
        'irate' : irate,
        'esyn_tau' : esyn_tau,
        'isyn_tau' : isyn_tau,
        'e_interval' : e_interval,
        'i_interval' : i_interval
}

par = {
        'gbar_naf_somatic': 9.0,
        'gbar_naf_axonal': 9.0
       }

       
nmda = {
        'gNMDAmax': gNMDAmax ,
        'gNMDAmax_plateau': gNMDAmax_plateau,
        'tcon': tNMDAon,
        'tcoff': tNMDAoff,
        'ratio': ratio_glutamate_syn
       }       
       
exp2syn = {
        'e': e_esyn ,
        'tau1': tau1_exp2syn ,
        'tau2': tau2_exp2syn
       }

#----------------------------------------------#
#      6. What to record       #
#----------------------------------------------#

#species_to_plot = ['EGF', 'EGFR', 'L_EGFR', 'SHC_p', 'SHC_phospho_cplx',
#                   'SHC_p_Sos_Grb2', 'GDP_Ras', 'Sos_Ras_GEF_cplx',
#                   'GTP_Ras', 'Raf_p_GTP_Ras', 'Raf_p_GTP_Ras_1_cplx',
#                   'MAPKK', 'MAPKK_p', 'MAPK', 'MAPK_p', 'MAPK_p_p',
#                   'MAPK_p_p_cplx', 'MAPK_p_p_feedback_cplx']            

species_to_plot = ['EGF']#, 'MAPK_p', 'MAPK_p_p']
pERK = ['MAPK_p', 'MAPK_p_p']

cascade_species = [
'EGF',
'PKC_Ca', 
'PKC_DAG_AA_p', 
'PKC_Ca_AA_p',  
'PKC_Ca_memb_p', 
'PKC_DAG_memb_p', 
'PKC_basal_p', 
'PKC_AA_p', 
'PKC_Ca_DAG', 
'PKC_DAG', 
'PKC_DAG_AA',
'PKC_cytosolic', 
'PLA2_cytosolic', 
'PLA2_Ca_p',
'PIP2_PLA2_p', 
'PIP2_Ca_PLA2_p', 
'DAG_Ca_PLA2_p', 
'PLA2_p_Ca', 
'PLA2_p', 
'Arachidonic_Acid', 
'PLC', 
'PLC_Ca', 
'PLC_Ca_Gq', 
'PLC_Gq',
'DAG', 
'IP3',
'MAPK_p_p', 
'craf_1',
'craf_1_p',
'MAPKK', 
'MAPK', 
'craf_1_p_p', 
'MAPK_p', 
'MAPKK_p_p', 
'MAPKK_p', 
'Raf_p_GTP_Ras', 
'craf_1_p_ser259', 
'inact_GEF', 
'GEF_p',
'GTP_Ras', 
'GDP_Ras',
'GAP_p',
'GAP', 
'inact_GEF_p', 
'CaM_GEF', 
'EGFR', 
'L_EGFR', 
'Internal_L_EGFR',
'SHC_p_Sos_Grb2', 
'SHC', 
'SHC_p', 
'Sos_p_Grb2',
'Grb2', 
'Sos_Grb2',
'Sos_p', 
'Sos', 
'SHC_p_Grb2_clx', 
'PLC_g', 
'PLC_g_p',
'Ca_PLC_g',
'Ca_PLC_g_p', 
'PLCg_basal', 
'MKP_1', 
'PPhosphatase2A', 
'PKC_act_raf_cplx', 
'PKC_inact_GAP_cplx', 
'PKC_act_GEF_cplx', 
'kenz_cplx', 
'kenz_cplx_1', 
'kenz_cplx_2', 
'kenz_cplx_3', 
'kenz_cplx_4', 
'PLC_Ca_cplx', 
'PLCb_Ca_Gq_cplx', 
'MAPK_p_p_cplx', 
'MAPK_p_p_feedback_cplx', 
'phosph_Sos_cplx', 
'MAPKKtyr_cplx', 
'MAPKKthr_cplx', 
'Raf_p_GTP_Ras_1_cplx', 
'Raf_p_GTP_Ras_2_cplx', 
'basal_GEF_activity_cplx', 
'GEF_p_act_Ras_cplx', 
'GAP_inact_Ras_cplx', 
'CaM_GEF_act_Ras_cplx', 
'Ca_PLC_g_phospho_cplx',
'SHC_phospho_cplx', 
'Sos_Ras_GEF_cplx', 
'PLC_g_phospho_cplx', 
'MKP1_tyr_deph_cplx', 
'MKP1_thr_deph_cplx', 
'craf_dephospho_cplx', 
'MAPKK_dephospho_cplx', 
'MAPKK_dephospho_ser_cplx', 
'craf_p_p_dephospho_cplx', 
'deph_raf_ser259_cplx']

#species_to_plot = cascade_species