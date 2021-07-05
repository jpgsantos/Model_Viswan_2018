# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 16:36:50 2017

@author: daniel
"""

from neuron import h
import experiment as e
import parameters as p
import time
import random as rnd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import json

class EGF_Experiment(e.Experiment):
    
    def __init__(self, exptype, cell):
        super(EGF_Experiment, self).__init__()
        self.exptype = exptype
        self.cell = cell
            
        if (str(type(cell))).find('MSN') != -1:
            self.celltype = 'MSN'
        
    def insert_synapses(self, syntype, syn_loc = [], deterministic = 0, 
                        num_syns = p.distributed_input_size, add_spine = 0, on_spine = 0):
                
        if syntype == 'input_syn':
            syn_loc = []            
            for i in range(0, num_syns):
                syn_loc.append([rnd.randint(0, len(self.cell.dendlist)-1), rnd.uniform(0,1)])
            
            if deterministic == 1:
                spike_time = []
                for i in range(0, len(syn_loc)):
                    spike_time.append(rnd.uniform(p.distributed_input_start, p.distributed_input_end))
            
            counter = 0
            for loc in syn_loc:
                counter += 1
                syn = self.cell.insert_synapse('glutamate', self.cell.dendlist[loc[0]], loc[1], add_spine = add_spine, on_spine = on_spine)
                if deterministic == 1:                
                    self.add_input_generator(syn, syntype, deterministic = deterministic, numsyn = counter, start = spike_time[counter-1])
                else:
                    self.add_input_generator(syn, syntype, deterministic = deterministic, numsyn = counter)
        
        elif syntype == 'MSN':
            for dend in self.cell.dendlist:
                # Synapses according to Cheng et al. Experimental Neurobiology, 147:287-298 (1997)                
                unit_length = 20.0 # Values reported in units per 20 microns in the study
                [exc_mean, exc_sem, inh_mean, inh_sem] = self.synapse_distribution(self.celltype, dend)
                # Insert excitatory synapses in this section, 
                # This contains both an exponential and an NMDA synapse.                
                syntype = 'glutamate'
                if dend.nseg <= exc_mean:
                    freq_multiplier = exc_mean/dend.nseg
                    step = 1/dend.nseg
                    for i in range(0, dend.nseg):
                        pos = (i + (i+1))*step/2
                        self.helper_insert(syntype, pos, dend, freq_multiplier)
                else:                        
                    num_exc_syn = int(dend.L/unit_length * rnd.gauss(exc_mean,exc_sem))
                    freq_multiplier = 1.0                     
                    for i in range(0,num_exc_syn):
                        pos = rnd.uniform(0,1)
                        self.helper_insert(syntype, pos, dend, freq_multiplier)                        
                
                # Insert inhibitory synapses in this section                
                syntype = 'inhexpsyn'
                if dend.nseg <= inh_mean:
                    freq_multiplier = inh_mean/dend.nseg
                    step = 1/dend.nseg
                    for i in range(0, dend.nseg):
                        pos = (i + (i+1))*step/2
                        self.helper_insert(syntype, pos, dend, freq_multiplier)
                else:        
                    num_inh_syn = int(dend.L/unit_length * rnd.gauss(inh_mean,inh_sem))                     
                    freq_multiplier = 1.0                    
                    for i in range(0,num_inh_syn):
                        pos = rnd.uniform(0,1)
                        self.helper_insert(syntype, pos, dend, freq_multiplier)  
        
        elif syntype in ['plateau_cluster']:
            if syntype == 'plateau_cluster':
                syntype = 'glutamate_ica_nmda'
            if deterministic == 0:
                if (type(syn_loc) == list):
                    for list_element in syn_loc:
                        if type(list_element) == list:
                            llim = list_element[1]; ulim = list_element[2];
                            loc = list_element[0]                        
                        else: # Now type(list_element) == int:
                            llim = 0.85; ulim = 1.0;
                            loc = list_element
                        for i in range(0, num_syns):
                            pos = rnd.uniform(llim, ulim)
                            syn = self.cell.insert_synapse(syntype, self.cell.dendlist[loc], pos, 
                                                           add_spine = add_spine, on_spine = on_spine)
                            self.add_input_generator(syn, syntype)
                            print(self.cell.dendlist[loc].name(), "%f" % (h.distance(llim, sec = self.cell.dendlist[loc])))
                
                elif (type(syn_loc) == dict):
                    for ind, loc in enumerate(syn_loc['loc']):
                        for i in range(0, num_syns):
                            syn = self.cell.insert_synapse(syntype, self.cell.dendlist[loc], 
                                                           syn_loc['pos'][ind], add_spine = add_spine, on_spine = on_spine)
                            self.add_input_generator(syn, syntype, 1.0 , start = syn_loc['start'][ind], 
                                                     end = syn_loc['end'][ind])
    
            elif deterministic == 1:
                if (type(syn_loc) == list):
                    for list_element in syn_loc: 
                        syn_step = 1.0/num_syns
                        for i in range(0, num_syns):
                            pos = p.cluster_end_pos - (p.cluster_end_pos - p.cluster_start_pos)*i*syn_step
                            syn = self.cell.insert_synapse(syntype, self.cell.dendlist[list_element], 
                                                           pos, add_spine = add_spine, on_spine = on_spine)
                            self.add_input_generator(syn, syntype, deterministic = deterministic, numsyn = i)
                
    def add_input_generator(self, syn, syntype, freq_multiplier = 1, 
                            start = p.plateau_burst_start, end = p.plateau_burst_end, deterministic = 0, numsyn = 1):
        
        if deterministic == 1:
            noise = 0
            start = start + numsyn*p.deterministic_interval
            number = p.num_deterministic_spikes
            interval = p.deterministic_interval
            weight = p.gAMPAmax_plateau
            
            if syntype in ['input_syn']:
                start = start
                weight = p.g_expsyn_max
            elif syntype in ['inhexpsyn_plateau']:
                weight = p.gGABAmax_plateau
                start = start + numsyn*2
                number = 3
                end = end    
            
        elif deterministic == 0:
            noise = 1
            if syntype in ['inhexpsyn']:
                start = 0;
                end = p.simtime
                number = (end-start) * p.irate * freq_multiplier
                interval = p.i_interval/freq_multiplier
                weight = p.g_inhexpsyn_max
                
            elif syntype in ['inhexpsyn_plateau']:
                start = p.inhibitory_burst_start
                end = p.inhibitory_burst_end
                number = (end-start) * p.inhibitory_syn_rate
                interval = p.inhibitory_syn_interval
                weight = p.gGABAmax_plateau

            elif syntype in ['expsyn_plateau', 'plateau_cluster', 
            'glutamate_ica_nmda']:
                number = (end-start) * p.plateau_syn_rate
                interval = p.plateau_syn_interval              
                weight = 1.0
                if syntype in ['expsyn_plateau']:
                    weight = p.gAMPAmax_plateau
                elif syntype in ['plateau_cluster']:
                    weight = p.gNMDAmax_plateau
                    
            elif syntype in ['input_syn']:
                start = p.distributed_input_start
                end = p.distributed_input_end
                number = (end-start) * p.distributed_input_rate
                interval = p.distributed_input_interval
                weight = p.gAMPAmax
                
            elif syntype in ['expsyn','glutamate']:
                start = 0;
                end = p.simtime
                number = (end-start) * p.erate * freq_multiplier
                interval = p.e_interval/freq_multiplier
                weight = p.g_expsyn_max
                
        gen = h.NetStim(0.5, sec = self.presyn)
        gen.seed(int(time.time() + rnd.randint(1,10**7)))
        gen.start = start
        gen.noise = noise
        gen.number = number        
        gen.interval = interval
        
        nc = h.NetCon(gen, syn.obj)
        nc.delay = 0
        nc.weight[0] = weight
        
        if syntype in ['inhexpsyn']:
            self.istim.append(gen) 
            self.inc.append(nc)      

        elif syntype in ['inhexpsyn_plateau']:
            self.istim.append(gen) 
            self.inc.append(nc)      

        elif syntype in ['expsyn', 'plateau_cluster', 'input_syn', 
        'expsyn_plateau', 'glutamate', 'glutamate_ica_nmda']:
            self.estim.append(gen)
            self.enc.append(nc)
                 
            nc = h.NetCon(gen, syn.obj)
            nc.delay = 0
            nc.weight[0] = weight
            self.estim.append(gen)
            self.enc.append(nc)
            
    def set_up_recording(self, dend_record_list = [], record_step = p.record_step):
        self.dend_record_list = dend_record_list        
        self.vdlist = []
        self.tout = h.Vector()
        self.tout.record(h._ref_t, record_step)
        self.tv = h.Vector()
        self.tv.record(h._ref_t, p.v_record_step)
        self.cali = []
        self.cali_dend = []
        self.cai_nmda = []
        self.cai = []            
        self.vspine = []

        self.vs = h.Vector()
        self.vs.record(self.cell.somalist[0](0.5)._ref_v, p.v_record_step)
                                      
        if self.exptype == 'single_sec':
            self.species = []            
            indices = []
            for s in p.cascade_species:
                indices.append(p.cascade_species.index(s))
            
            for i in range(0,len(p.cascade_species)):
                self.species.append(h.Vector())
                if not self.cell.opt_flag:
                    cmd = 'self.cell.somalist[0](0.5)._ref_' + p.cascade_species[indices[i]] + '_Viswan_2018'                                                
                else:
                    cmd = 'self.cell.somalist[0](0.5)._ref_' + p.cascade_species[indices[i]] + '_Viswan_2018_opt'                                                
                self.species[-1].record(eval(cmd), record_step)
            
        if self.exptype == 'record_ca' or self.exptype == 'record_ca_opt':
            self.species = []
            self.pERK_species = []
            self.integral = []
            indices = []
            indices_pERK = []
            self.pERKratio = h.Vector()
            if self.exptype == 'record_ca' and p.pERK != []:
                self.pERKratio.record(self.cell.somalist[0](0.5)._ref_pERK1_2_ratio1_Viswan_2018, record_step)
            elif self.exptype == 'record_ca_opt' and p.pERK != []:
                self.pERKratio.record(self.cell.somalist[0](0.5)._ref_pERK1_2_ratio1_Viswan_2018_opt, record_step)
            for s in p.species_to_plot:
                indices.append(p.cascade_species.index(s))

            for s in p.pERK:
                indices_pERK.append(p.cascade_species.index(s))

            for d in dend_record_list:
                self.vdlist.append(h.Vector())
                self.vdlist[-1].record(self.cell.dendlist[d](p.pos)._ref_v, record_step)
                sec = self.cell.dendlist[d]
#                record_spinelist = [s for s in self.cell.spines if s.parent == self.cell.dendlist[d]]                                        
#                spine = record_spinelist[-1]                    
                self.cai.append(h.Vector())
                self.cali.append(h.Vector())
                self.cai_nmda.append(h.Vector())
                self.cali_dend.append(h.Vector())
#                self.vspine[-1].record(spine.head(0.5)._ref_v, record_step)                    
                self.cai[-1].record(sec(p.pos)._ref_cai, record_step)
                self.cali[-1].record(sec(p.pos)._ref_cali, record_step)
                self.cai_nmda[-1].record(sec(p.pos)._ref_ca_nmdai, record_step)
                self.cali_dend[-1].record(sec(p.pos)._ref_cali, record_step)
            
            for i in range(0,len(p.species_to_plot)):
                self.species.append(h.Vector())
                if self.exptype == 'record_ca':
                    cmd = 'self.cell.somalist[0](0.5)._ref_' + p.cascade_species[indices[i]] + '_Viswan_2018'                                                
                elif self.exptype == 'record_ca_opt':
                    cmd = 'self.cell.somalist[0](0.5)._ref_' + p.cascade_species[indices[i]] + '_Viswan_2018_opt'                                                
                self.species[-1].record(eval(cmd), record_step)

            for i in range(0,len(p.pERK)):
                self.pERK_species.append(h.Vector())
                if self.exptype == 'record_ca':
                    cmd = 'self.cell.somalist[0](0.5)._ref_' + p.cascade_species[indices_pERK[i]] + '_Viswan_2018'                                                
                elif self.exptype == 'record_ca_opt':
                    cmd = 'self.cell.somalist[0](0.5)._ref_' + p.cascade_species[indices_pERK[i]] + '_Viswan_2018_opt'                                                
                self.pERK_species[-1].record(eval(cmd), record_step)

                            
            self.w_ampa = []
            
    def plot_results(self):
#        sns.set(font_scale = 2.0)
        sns.set_style("whitegrid")                
        plt.style.use('ggplot')
                   
        if self.exptype == 'single_sec':

            figs = []    
            axes = []
#            for i in range(0,len(p.species_to_plot)):
#                figs.append(plt.figure())
#                axes.append(figs[-1].add_subplot(111))
#                axes[-1].set_ylabel(p.species_to_plot[i] + ' (uM)')
#                axes[-1].set_xlabel('t')
#                axes[-1].plot(self.tout, self.species[i])
            plt.show()

        if self.exptype == 'record_ca' or self.exptype == 'record_ca_opt':
            fig_vs = plt.figure()
            ax_vs = fig_vs.add_subplot(111)
            ax_vs.set_ylabel('Vs')
            ax_vs.set_xlabel('t')
            ax_vs.plot(self.tv, self.vs)
            
            ss = np.asarray(self.pERK_species[0].to_python()) \
               + np.asarray(self.pERK_species[1].to_python()) 

            fig_pERK = plt.figure()
            ax_pERK = fig_pERK.add_subplot(111)
            ax_pERK.set_ylabel('pERK1_2_ratio1')
            ax_pERK.set_xlabel('t')            
            ax_pERK.plot(self.tout, ss)

            fig_pERKratio = plt.figure()
            ax_pERKratio = fig_pERKratio.add_subplot(111)
            ax_pERKratio.set_ylabel('pERKratio')
            ax_pERKratio.set_xlabel('t')            
            ax_pERKratio.plot(self.tout, self.pERKratio)

            
#           1 ax_vd.plot(self.tout, self.vdlist[i])
                    
#            legend_list = []
#            for i in self.dend_record_list:
#                legend_list.append('dend[%d](%.2f) = %.2f um' % (i, p.pos, h.distance(p.pos, sec = self.cell.dendlist[i]) ))
                
#            fig_cai_nmda = plt.figure(); 
#            ax_cai_nmda = fig_cai_nmda.add_subplot(111);
#            ax_cai_nmda.set_ylabel('[Ca]_NMDA')
#            ax_cai_nmda.set_xlabel('t')
#            plt.hold(True)
#            for i in range(0,len(self.cai_nmda)):
#                ax_cai_nmda.plot(self.tout, self.cai_nmda[i])

#            if self.cell.spines != []:
#                fig_vspine = plt.figure()
#                ax_vspine = fig_vspine.add_subplot(111)
#                ax_vspine.set_ylabel('Vspine')
#                ax_vspine.set_xlabel('t')
#                for v in self.vspine:
#                    ax_vspine.plot(self.tout, v)
#            
#            legend_list = []
#            for i in self.dend_record_list:
#                legend_list.append('dend[%d](%.2f) = %.2f um' % (i, p.pos,  h.distance(p.pos, sec = self.cell.dendlist[i]) ))
            
#            ax_vd.legend(legend_list);  
#            ax_cai_nmda.legend(legend_list);

            figs = []    
            axes = []
            for i in range(0,len(p.species_to_plot)):
                figs.append(plt.figure())
                axes.append(figs[-1].add_subplot(111))
                axes[-1].set_ylabel(p.species_to_plot[i] + ' (uM)')
                axes[-1].set_xlabel('t')
                axes[-1].plot(self.tout, self.species[i])
            
            plt.show()
            return fig_vs

    def plot_helper(self, plot_what):
        self.synlist = self.get_synapse_list('adaptive_glutamate2')
        self.synlist.sort(key = lambda f: f.sec.name())

        xlabel = 't'
        if plot_what == 'wampa':
            yval = [s.ref_var_ampa for s in self.synlist]
            ylabel = 'wampa'
            yliml = p.LTD_factor*p.gAMPAmax_plateau
            ylimu = p.LTP_factor*p.gAMPAmax_plateau
        elif plot_what == 'wnmda':
            yval = [s.ref_var_nmda for s in self.synlist]
            ylabel = 'wnmda'
            yliml = p.LTD_factor*p.gNMDAmax_plateau
            ylimu = p.LTP_factor*p.gNMDAmax_plateau
        elif plot_what == 'ca_nmda':
            yval = self.cai_nmda_in_syns
            ylabel = 'ca_nmda'
        elif plot_what == 'cali':
            yval = self.cali_in_syns
            ylabel = 'cali'
        
        fig = []; axes = []; legend = []
        for i in range(0,len(self.synlist)):
            if i==0 or (self.synlist[i].sec.name() != self.synlist[i-1].sec.name()):                   
                fig.append(plt.figure())
                axes.append(fig[-1].add_subplot(111))
                plt.hold(True)
                axes[-1].set_xlabel(xlabel) 
                axes[-1].set_ylabel(ylabel)
                if plot_what == 'wampa' or plot_what == 'wnmda':
                    axes[-1].set_ylim(yliml, ylimu)
                    
                if not (i == 0):
                    axes[-2].legend(legend)
                    legend = []
            
            axes[-1].plot(self.tout, yval[i])
            
            string = '%s(%.2f) = %.2f um' % (self.synlist[i].sec.name(), self.synlist[i].pos,
                    h.distance(self.synlist[i].pos, sec = self.synlist[i].sec) )
            legend.append(string)
            
        return fig, axes
           
    def simulate(self, simtime = p.simtime, parallel = False):
        start = time.time()
        gmtime = time.gmtime(start)
        print("Starting simulation... %d:%d:%d" %(gmtime.tm_hour + 1, gmtime.tm_min, gmtime.tm_sec))
        
        if not parallel:        
            h.load_file("stdrun.hoc")
            h.init()
            h.tstop = simtime        
            p.simtime = simtime       
            fih1 = h.FInitializeHandler((self.seti_print_status))
            fih2 = h.FInitializeHandler((self.seti_steady_states))
#            fih3 = h.FInitializeHandler((self.seti_EGF, p.EGF_input_start))
#            cvode = h.CVode()
#            cvode.active(1)
#
#            update_points = np.arange(0, p.simtime + p.simtime/100., p.simtime/100. )        
#            while h.t < simtime:
#                if h.t > update_points[0]:
#                    print(update_points[0])
#                    update_points = update_points[1:len(update_points)]
#                h.fadvance()
            h.run()

        end = time.time()
        print("It took %.2f hours or %.4f seconds to simulate." % ((end-start)/3600 , (end-start)))

    def seti_steady_states(self):
        if self.exptype == 'record_ca':
            filename = 'steady_states.dat'
        else:
            filename = 'steady_states_opt.dat'
        with open(filename, 'r', encoding = 'utf-8') as f:
            to_read = json.load(f)
            result = json.loads(to_read)
            steady_states = result['steady_states']    
            for i in range(0,len(steady_states)):
                if self.exptype == 'record_ca':
                    cmd = 'self.cell.somalist[0](0.5).' + p.cascade_species[i] + '_Viswan_2018' + '=' + str(steady_states[i])
                    cmd_print = 'self.cell.somalist[0](0.5).' + p.cascade_species[i] + '_Viswan_2018'                
                elif self.exptype == 'record_ca_opt':
                    cmd = 'self.cell.somalist[0](0.5).' + p.cascade_species[i] + '_Viswan_2018_opt' + '=' + str(steady_states[i])
                    cmd_print = 'self.cell.somalist[0](0.5).' + p.cascade_species[i] + '_Viswan_2018_opt'                
                exec(cmd)
                print(p.cascade_species[i], '=', eval(cmd_print) )                
                    
    def seti_print_status(self):
        update_points = np.arange(0, p.simtime, p.simtime/100. )
        for t in update_points:
            h.cvode.event(t, self.print_status)
    
    def print_status(self):
        print("At time t %f, total simtime %f." % (h.t, p.simtime))

    def seti_EGF(self, time):
        h.cvode.event(time, (self.set_EGF_input, p.EGF_conc))
    
    def set_EGF_input(self, value):
        h.EGF = value
        if not value == 0.0:
            h.cvode.event(h.t + p.EGF_input_length, (self.set_EGF_input, 0.0))
        
    def calculate_aoc(self, vector):
        v = np.array(vector.to_python())
        return sum(v) - v[0]*(len(v))
        
    def get_synapse_list(self, syntype):
        synlist = []
        for syn in self.cell.esyn:
            if syn.type == syntype:
                synlist.append(syn)
        
        return synlist
        
    def insert_IClamp(self, sec, pos, num_clamps):
        self.iclamp = []
        for i in range(0,num_clamps):                   
            self.iclamp.append(h.IClamp(pos, sec=sec))
            self.iclamp[-1].amp = p.iclamp_amp 
            self.iclamp[-1].delay = p.plateau_burst_start + p.iclamp_delay + i*p.iclamp_periodic_delay
            self.iclamp[-1].dur = p.iclamp_dur