
import GillespieSim as gs
import numpy as np
import random


                                # _____ Syntaxis _____ #

# species =       {'species_name': quantity_of_molecules}
# propensities =  {'propensities_name': propensity_value}
# reactions =     {'reaction_name': ['propensity', 'species implyied']}
# transcription = {'reaction_name':{'create': ['species'] , 'destroy': [None]}}

                                # _____         _____ #
    
central_dogma = gs.ReactionModel(
                                    species= {
                                            'dna': 1, 
                                            'rna': 0, 
                                            'protein': 0, 
                                            'complex': 0
                                            },
                                    propensities={
                                                'trc_c':    0.3, 
                                                'trl_c':    10.02, 
                                                'pdeg_d':   np.log(2)/30,
                                                'rdeg_d':   np.log(2)/2.5, 
                                                'cmplx_c':  0.0006,
                                                'cmplx_d':  0.3
                                            },
                                    reactions={
                                                'transcription':   ['trc_c', 'dna'], 
                                                'translation':     ['trl_c', 'rna'], 
                                                'p_degradation':   ['pdeg_d', 'protein'], 
                                                'rna_degradation': ['rdeg_d', 'rna'],
                                                'complex_c':       ['cmplx_c', 'dna', 'protein'],
                                                'complex_d':       ['cmplx_d', 'complex']
                                            },
                                    q={
                                        'transcription':   {'create': ['rna'],            'destroy': [None]}, 
                                        'translation':     {'create': ['protein'],        'destroy': [None]}, 
                                        'p_degradation':   {'create': [None],             'destroy': ['protein']}, 
                                        'rna_degradation': {'create': [None],             'destroy': ['rna']},
                                        'complex_c':       {'create': ['complex'],        'destroy': ['dna', 'protein']},
                                        'complex_d':       {'create': ['dna', 'protein'], 'destroy': ['complex']}
                                    } 
                                    
                                )  

birth_death = gs.ReactionModel(
    
                                species= {
                                            'protein': 0
                                            },
                                propensities= {
                                                'kr': 100, 
                                                'gamma':10
                                                },

                                reactions= {
                                            'translation':['kr'], 
                                            'degradation':['gamma', 'protein']
                                            },
                                q= {
                                    'translation': {'create': ['protein'],  'destroy': [None]}, 
                                    'degradation': {'create': [None],       'destroy': ['protein']}
                                    }
                                    )

lamp_1 = gs.ReactionModel(
                            species = {
                                        'dna': 5,
                                        'amplicons': 0,
                                        'primers': 10000
                                        },

                            propensities = {
                                            'pc': 0.00007
                                            },

                            reactions = {
                                        'amplification': ['pc', 'dna', 'primers']
                                        },

                            q = {
                                'amplification': {'create': ['dna', 'amplicons'], 'destroy': ['primers']}
                                }
                        )

lamp_2 = gs.ReactionModel(
                            species = {
                                            'in_primer_a': random.randint(900, 1000),
                                            'ex_primer_a': random.randint(120, 150),
                                            'gene': random.randint(3, 5),
                                            'amplicon': 0},

                            propensities = {
                                                'gene_c': 0.00039,
                                                'amplicon_c': 0.000185,
                                                'amplicon_amp': 0.0025},

                            reactions = {
                                            'amplification':    ['gene_c',          'in_primer_a',  'gene'],
                                            'loop_amp':         ['amplicon_c',      'ex_primer_a',  'gene'],
                                            'amp_amp':          ['amplicon_amp',    'ex_primer_a',  'amplicon']},

                            q = {
                                    'amplification':    {'create': ['gene'],     'destroy': ['in_primer_a']},
                                    'loop_amp':         {'create': ['amplicon'], 'destroy': ['ex_primer_a', 'gene']},
                                    'amp_amp':          {'create': ['amplicon'], 'destroy': ['ex_primer_a']}}
                        )

lamp_3 = gs.ReactionModel(
                            species = {
                                            'f_inner_primer': 4140,
                                            'b_inner_primer': 4140,
                                            'f_outer_primer': 1035,
                                            'b_outer_primer': 1035,
                                            'i_gene': 3,
                                            'f_gene':0,
                                            'r_gene': 0,
                                            'amplicon': 0
                                            },
                                            
                            propensities = {
                                                'krg': 0.00035,
                                                'kcamp': 0.0004
                                                },
                            reactions = {
                                            'gene_amp': ['krg', 'f_inner_primer', 'b_inner_primer','i_gene'],
                                            'loop_amplicon': ['kcamp', 'f_outer_primer','b_outer_primer', 'r_gene']
                                            },
                            q = {
                                    'gene_amp': {'create': ['r_gene'],  'destroy': ['inner_primer']}, 
                                    'loop_amplicon': {'create': ['amplicon'],  'destroy': ['outer_primer']}
                            })

repressilator = gs.ReactionModel(
                                species = {
                                                'tetR': 550,
                                                'lacI': 1600,
                                                'alphacl': 1000,
                                                'GFP': 100
                                                },

                                propensities = {
                                                'kd_tetR': 0.02,
                                                'kd_lacI': 0.02,
                                                'kd_alphacl': 0.02,
                                                'kc_GFP': 0.05,
                                                'kd_GFP': 1/100000
                                                },
                                
                                reactions = {
                                                'PdtetR':       ['kd_tetR', 'alphacl'],
                                                'PdlacI':       ['kd_lacI', 'tetR'],
                                                'Pdalphacl':    ['kd_alphacl', 'lacI'],
                                                'GFP_cre' :     ['kc_GFP', 'lacI'],
                                                'GFP_des':      ['kd_GFP', 'tetR', 'GFP']
                                                },
                                
                                q = {
                                        'PdtetR':       {'create': ['lacI'],    'destroy': ['tetR']},
                                        'PdlacI':       {'create': ['alphacl'], 'destroy': ['lacI']},
                                        'Pdalphacl':    {'create': ['tetR'],    'destroy': ['alphacl']},
                                        'GFP_cre':      {'create': ['GFP'],     'destroy': [None]},
                                        'GFP_des':      {'create': [None],      'destroy': ['GFP']}
                                        }
                                )

repressilator = gs.ReactionModel(
                                species = {
                                                'tetR': 0,
                                                'lacI': 0,
                                                'alphacl': 0,
                                                'GFP': 0,
                                                'complex_tetR': 0,
                                                'complex_lacI': 0,
                                                'complex_alphacl': 0,
                                                'complex_GFP': 0,
                                                'promoter_tetR': 1,
                                                'promoter_lacI': 1,
                                                'promoter_alphacl': 1,
                                                'promoter_GFP': 1
                                                },

                                propensities = {
                                                'kc_tetR': 10.02,
                                                'kc_lacI': 10.02,
                                                'kc_alphacl': 10.02,
                                                'kd_tetR': np.log(2)/30,
                                                'kd_lacI': np.log(2)/30,
                                                'kd_alphacl': np.log(2)/30,
                                                'kc_GFP': 10.02,
                                                'kd_GFP': np.log(2)/60,
                                                'ka_c_tetR': 0.00001,
                                                'ka_c_lacI': 0.00001,
                                                'ka_c_alphacl': 0.00001,
                                                'ka_c_GFP': 0.00001,
                                                'kd_c_tetR': 0.005,
                                                'kd_c_lacI': 0.005,
                                                'kd_c_alphacl': 0.005,
                                                'kd_c_GFP': 0.005
                                                },


                                reactions = {
                                                'PctetR':       ['kc_tetR', 'promoter_tetR'],
                                                'PclacI':       ['kc_lacI', 'promoter_lacI'],
                                                'Pcalphacl':    ['kc_alphacl', 'promoter_alphacl'],
                                                'PdtetR':       ['kd_tetR', 'tetR'],
                                                'PdlacI':       ['kd_lacI', 'lacI'],
                                                'Pdalphacl':    ['kd_alphacl', 'alphacl'],
                                                'GFP_cre' :     ['kc_GFP', 'promoter_GFP'],
                                                'GFP_des':      ['kd_GFP', 'GFP'],
                                                'Ac_tetR':      ['ka_c_tetR', 'lacI'],
                                                'Ac_lacI':      ['ka_c_lacI', 'alphacl'],
                                                'Ac_alphacl':   ['ka_c_alphacl', 'tetR'],
                                                'Ac_GFP':       ['ka_c_GFP', 'tetR'],
                                                'Dc_tetR':      ['kd_c_tetR'],
                                                'Dc_lacI':      ['kd_c_lacI'],
                                                'Dc_alphacl':   ['kd_c_alphacl'],
                                                'Dc_GFP' :      ['kd_c_GFP']},
                                
                                q = {   'PctetR':       {'create': ['tetR'],    'destroy': [None]},
                                        'PclacI':       {'create': ['lacI'],    'destroy': [None]},
                                        'Pcalphacl':    {'create': ['alphacl'], 'destroy': [None]},
                                        'GFP_cre':      {'create': ['GFP'],     'destroy': [None]},
                                        'PdtetR':       {'create': [None],      'destroy': ['tetR']},
                                        'PdlacI':       {'create': [None],      'destroy': ['lacI']},
                                        'Pdalphacl':    {'create': [None],      'destroy': ['alphacl']},           
                                        'GFP_des':      {'create': [None],      'destroy': ['GFP']},
                                        'Ac_tetR' :     {'create': ['complex_tetR'],            'destroy': ['lacI', 'promoter_tetR']},
                                        'Ac_lacI' :     {'create': ['complex_lacI'],            'destroy': ['alphacl', 'promoter_lacI']},
                                        'Ac_alphacl' :  {'create': ['complex_alphacl'],         'destroy': ['tetR', 'promoter_alphacl']},
                                        'Ac_GFP':       {'create': ['complex_GFP'],             'destroy': ['tetR', 'promoter_GFP']},
                                        'Dc_tetR' :     {'create': ['lacI', 'promoter_tetR'],   'destroy': ['complex_tetR']},
                                        'Dc_lacI' :     {'create': ['alphacl', 'promoter_lacI'],'destroy': ['complex_lacI']},
                                        'Dc_alphacl' :  {'create': ['tetR', 'promoter_alphacl'],'destroy': ['complex_alphacl']},
                                        'Dc_GFP' :      {'create': ['tetR', 'promoter_GFP'],    'destroy': ['complex_GFP']}
                                        }
                                )


lacI_model = gs.ReactionModel(

        species = {
                'mRNA': 10,
                'Complex': 0,
                'DNA': 1,
                'Beta': 0,              
                'Intracellular_tmg': 0,
                'LacI_Tetramers': 0,   
                'Active_LacI': 0,      
                'Permease': 0,        
                'LacI_monomer': 0,
                'TMG': 25},  
        
        propensities = {
                        'kc_mRNA': 2,
                        'kd_mRNA': np.log(2)/2.5,
                        'kp_on': 1.5 * (6.666/666.666),
                        'kp_off': 1.5 * (0.04/666.666)},

        reactions = {
                'Transcription': ['kc_mRNA', 'DNA'],
                'mRNA_Degradation': ['kd_mRNA', 'mRNA'],
                'Complex_Creation': ['kp_on', 'DNA'],
                'Complex_Degradation' : ['kp_off']},
        
        q = {
        'Transcription':       {'create': ['mRNA'],       'destroy': [None]},
        'mRNA_Degradation':    {'create': [None],         'destroy': ['mRNA']},
        'Complex_Creation':    {'create': ['Complex'],    'destroy': ['DNA']},
        'Complex_Degradation': {'create': ['DNA'],        'destroy': ['Complex']}},
        
        math_model= {
                'Beta':              lambda tmg:                         0.00123 * pow(tmg, 0.6),
                'Intracellular_tmg': lambda beta, permease:              beta * permease,
                'LacI_Tetramers':    lambda monomer:                     monomer + (monomer/100),
                'Active_LacI':       lambda intracellular_tmg, tetramer: 1/(1 + (intracellular_tmg / 0.12)**2) * tetramer,
                'Permease':          lambda mrna, permease, tau:         (1000 * 0.00231 * mrna/1000) + (permease - (1000 * 0.00231 * mrna/1000)) * np.exp(-1000 * tau),
                'LacI_monomer':      lambda monomer, tau:                (15.4 - 0.0231 * monomer)*tau},

        interaction_map = {       
                        'Beta':                 ['TMG'],
                        'Intracellular_tmg':    ['Beta', 'Permease'],
                        'LacI_Tetramers':       ['LacI_monomer'],
                        'Active_LacI':          ['Intracellular_tmg', 'LacI_Tetramers'],
                        'Permease':             ['mRNA', 'Permease', 'Tau'],
                        'LacI_monomer':         ['LacI_monomer', 'Tau']}
                        )
