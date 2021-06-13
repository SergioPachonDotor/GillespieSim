import numpy as np 
import pandas as pd
import copy
import warnings
import csv
import os
from tqdm import tqdm
import random

class Gillespie:

    def __init__(self, tmax=10, sampling_time = 0.01,reaction_model = object, model_name='test', cells=0):
        self.sampling_time = sampling_time      # Sampling Time
        self.time = 0                           # Iitial Time
        self.reference_time = 0                 # Reference Time
        self.tmax = tmax                        # Total Time
        self.reaction_model = reaction_model    # Object of reaction model
        self._reaction_model_copy = copy.deepcopy(reaction_model)
        self.model_name = f'{model_name}.csv'
        self.tau = 0
        self.cells = cells

  
    def simulate(self):
        """ 
            Performs Cesar Nieto et-al Gillespie algorithm
            for chemical reactions simulations.
        """
        schema = self.create_schema()

        with open(self.model_name, mode='a', newline='') as f:

            writer = csv.DictWriter(f, fieldnames= schema, delimiter='|')

            for cell in tqdm(range(1, self.cells+1)):
                while self.time < self.tmax:

                    propensities = self.reaction_model.calculate_propensities()
                    propenisties_list = list(propensities.values())
                    possible_reactions = None

                    warnings.simplefilter("ignore")

                    reaction_times = lambda propensity, random_num: -(1/propensity) * np.log(random_num)
                    random_array = [np.random.rand() for k in range(len(propenisties_list))]

                    calculated_propensities = list(map(reaction_times, propenisties_list, random_array))
                    reactions_dict = dict(enumerate(list(self.reaction_model.show_reactions().keys())))
                    
                    self.tau = np.min(calculated_propensities)     # tau time
                    q = np.argmin(calculated_propensities)    # Reaction that occurs

                    if self.time + self.tau < self.reference_time:
                        #
                        # Integrate Functions
                        #
                        possible_reactions = self.reaction_model.show_q()[reactions_dict[q]]

                        create_species = possible_reactions['create']
                        [self.reaction_model.create(name=species) for species in create_species if species != None]

                        destroy_species = possible_reactions['destroy']
                        [self.reaction_model.destroy(name=species) for species in destroy_species if species != None]

                        self.time += self.tau
                        
                    else:
                        #
                        # Integrate Functions
                        #
                        self.time = self.reference_time  
                        model = self.model_data(cell=cell)
                        writer.writerow(model)

                        self.reference_time += self.sampling_time
                        model = None
                
                self.reaction_model = copy.deepcopy(self._reaction_model_copy)
                self.time = 0
                self.reference_time = 0
        f.close()

    def create_file(self):

        schema= self.create_schema()
        try: 
            os.remove(self.model_name)
        except:
            pass
        try:
            with open(self.model_name, mode='x') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=schema, delimiter='|')
                writer.writeheader()
                csvfile.close()
        except:
            print(f'File "{self.model_name}" already exist')
            pass

    def create_schema(self):
        schema = copy.deepcopy(self.reaction_model.show_species())
        schema = list(schema.keys())
        schema.append('time')
        schema.append('cell')
        return schema

    def model_data(self, cell):
        model = copy.deepcopy(self.reaction_model.show_species())
        model.update({'time': round(self.time, 2), 'cell': cell})
        return model

    def load_data(self):
        df = pd.read_csv(self.model_name, delimiter='|')
        return df

class ReactionModel:
    
    def __init__(self, reactions=dict,species=dict, propensities=dict, q=dict) -> None:
        self.reactions = reactions
        self.species = species 
        self.propensities = propensities
        self.q = q
    
    def show_species(self):
        return self.species

    def show_reactions(self):
        return self.reactions

    def show_propensities(self):
        return self.propensities
    
    def show_q(self):
        return self.q
    
    def calculate_propensities(self):
        
        complete_propensity = list(self.reactions.values())
        reaction_names = self.reactions.keys()
        species_copy = copy.copy(self.species)
        species_copy.update(self.propensities)

        species_calculus = []

        for i in range(len(complete_propensity)):
            propensities_array = []
            species_calculus.append(propensities_array)

            for k in complete_propensity[i]:
                propensities_array.append(species_copy[k])
                
        propensities_result = list(map(np.prod, species_calculus))
        
        return dict(zip(reaction_names, propensities_result))
    
    def create(self, name=str):
        self.species[name] += 1
    
    def destroy(self, name=str):
        self.species[name] -= 1
    
    def set_q(self):
        pass

class FunctionModel:
    def __init__(self) -> None:
        pass

if __name__ == '__main__':

    # _____ Syntaxis _____ #

    # species = {'species_name': quantity_of_molecules}
    # propensities = {'propensities_name': propensity_value}
    # reactions = {'reaction_name': ['propensity', 'species implyied']}
    # transcription = {'reaction_name':{'create': ['species'] , 'destroy': [None]}}

    # _____     _____ #

    species = {
                'dna': 1, 
                'rna': 0, 
                'protein': 0, 
                'complex': 0
                }

    propensities = {
                    'trc_c':    0.3, 
                    'trl_c':    10.02, 
                    'pdeg_d':   np.log(2)/30,
                    'rdeg_d':   np.log(2)/2.5, 
                    'cmplx_c':  0.0006,
                    'cmplx_d':  0.3
                }

    reactions = {
                    'transcription':   ['trc_c', 'dna'], 
                    'translation':     ['trl_c', 'rna'], 
                    'p_degradation':   ['pdeg_d', 'protein'], 
                    'rna_degradation': ['rdeg_d', 'rna'],
                    'complex_c':       ['cmplx_c', 'dna', 'protein'],
                    'complex_d':       ['cmplx_d', 'complex']
                }

    q = {
            'transcription':   {'create': ['rna'],            'destroy': [None]}, 
            'translation':     {'create': ['protein'],        'destroy': [None]}, 
            'p_degradation':   {'create': [None],             'destroy': ['protein']}, 
            'rna_degradation': {'create': [None],             'destroy': ['rna']},
            'complex_c':       {'create': ['complex'],        'destroy': ['dna', 'protein']},
            'complex_d':       {'create': ['dna', 'protein'], 'destroy': ['complex']}
        }  

    # model = ReactionModel(reactions=reactions, species=species, propensities=propensities, q=q)
    # Gillespie(tmax=720, sampling_time=0.1,reaction_model = model).simulate(); 

    species_2 = {
                    'protein': 0
                    }

    propensities_2 = {
                        'kr': 100, 
                        'gamma':10
                        }
    
    reactions_2 = {
                    'translation':['kr'], 
                    'degradation':['gamma', 'protein']
                    }
    
    q_2 = {
            'translation': {'create': ['protein'],  'destroy': [None]}, 
            'degradation': {'create': [None],       'destroy': ['protein']}
        }

    # model_2 = ReactionModel(reactions=reactions_2, species=species_2,propensities=propensities_2,q=q_2)
    # model_2_gillespie = Gillespie(tmax=10, sampling_time=0.01, reaction_model=model_2, model_name='birth_death')
    # model_2_gillespie.create_file()
    # model_2_gillespie.simulate()


    species_3 = {
        'dna': 5,
        'amplicons': 0,
        'primers': 10000,}

    propensities_3 = {
        'pc': 0.00007}

    reactions_3 = {
        'amplification': ['pc', 'dna', 'primers'],}

    q_3 = {
        'amplification': {'create': ['dna', 'amplicons'], 'destroy': ['primers']}}

    # model_3 = ReactionModel(reactions=reactions_3, species=species_3, propensities=propensities_3, q=q_3)
    # Gillespie(tmax=20, sampling_time=0.01, reaction_model=model_3).simulate()

    species_4 = {
                    'A': 10,
                    'B': 1,
                    'Z': 0,
                    'C': 0}

    propensities_4 = {
                        'kcza': 10,
                        'kcaa': 0.5,
                        'kcba': 0.5,
                        'kcbb': 10,
                        'kcc': 0.001,
                        'kdzb': 6,
                        'kdc': 0.5}

    reactions_4 = {
                    'cza': ['kcza', 'A'],
                    'cba': ['kcba', 'A'],
                    'caa': ['kcaa', 'A'],
                    'cbb': ['kcbb', 'B'],
                    'dab': ['kcc', 'A', 'B'],
                    'dzb': ['kdzb', 'B', 'Z'],
                    'dc' : ['kdc', 'C']}

    q_4 = {
            'cza': {'create' : ['Z'], 'destroy' : ['A']},
            'cba': {'create' : ['B'], 'destroy' : ['A']},
            'caa': {'create' : ['A'], 'destroy' : [None]},
            'cbb': {'create' : ['B'], 'destroy' : [None]},
            'dab': {'create' : ['C'], 'destroy' : ['A', 'B']},
            'dzb': {'create' : [None], 'destroy' : ['Z', 'B']},
            'dc' : {'create' : ['A', 'B'], 'destroy' : ['C']}}

    # model_4 = ReactionModel(
    #                         reactions=reactions_4, 
    #                         species=species_4, 
    #                         propensities=propensities_4, 
    #                         q=q_4)

    # Gillespie(tmax=60, sampling_time=0.8, reaction_model=model_4).simulate()

    species_5 = {
                    'in_primer_a': random.randint(900, 1000),
                    'ex_primer_a': random.randint(120, 150),
                    'gene': random.randint(3, 5),
                    'amplicon': 0}

    propensities_5 = {
                        'gene_c': 0.00039,
                        'amplicon_c': 0.000185,
                        'amplicon_amp': 0.0025}

    reactions_5 = {
                    'amplification':    ['gene_c',          'in_primer_a',  'gene'],
                    'loop_amp':         ['amplicon_c',      'ex_primer_a',  'gene'],
                    'amp_amp':          ['amplicon_amp',    'ex_primer_a',  'amplicon']}

    q_5 = {
            'amplification':    {'create': ['gene'],     'destroy': ['in_primer_a']},
            'loop_amp':         {'create': ['amplicon'], 'destroy': ['ex_primer_a', 'gene']},
            'amp_amp':          {'create': ['amplicon'], 'destroy': ['ex_primer_a']}}

    model_5 = ReactionModel(
                        reactions=reactions_5, 
                        species=species_5, 
                        propensities=propensities_5, 
                        q=q_5)

    sim = Gillespie(tmax=30, sampling_time=0.1, reaction_model=model_5, model_name='lamp_test', cells=1000)
    sim.create_file()
    sim.simulate()


    species_6 = {
                    'bst_polymerase': 1,
                    'dntp':  214120,
                    'in_primer_a': random.randint(900, 1000),
                    'ex_primer_a': random.randint(120, 150),
                    'gene': random.randint(3, 5),
                    'amplicon': 0}
                    
    propensities_6 = {
                        'amplicon_propensity': 39,
                        '': 0
                        }
    reactions_6 = {}
    q_5 = {}