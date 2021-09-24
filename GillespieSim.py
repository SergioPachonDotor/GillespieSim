import numpy as np
import pandas as pd
import copy
import warnings
import csv
import os
from tqdm import tqdm
from itertools import starmap

# ___________________________________________ 
# 
# Code Wrtitten By Sergio Andrés Pachón Dotor
# 2020
# e-mail: sap9827@gmail.com
# 
# ___________________________________________



class Gillespie:


    def __init__(
                    self, 
                    tmax=60, 
                    sampling_time = 0.1,
                    reaction_model = object, 
                    model_name='test', 
                    cells=1, 
                    time_stamp=False, 
                    division_time=18,
                    division_size=2,
                    birth_size=np.random.normal(loc=2/2, scale=0.1*(2/2)),
                    burst=15,
                    d_burst=1.5,
                    chunk_lenght = 100
                    ):

        # __ Parameters ___ #
        
        # Time Parameters
        self.sampling_time = sampling_time      # Sampling Time
        self.time = 0                           # Iitial Time
        self.reference_time = 0                 # Reference Time

        # Size Parameters
        self.mu = np.log(2)/division_time       # Mu parameter for Sizer Model
        self.division_size = division_size
        
        self.birth_size = birth_size
        self.size = 0

        # Cells
        self.cells = cells                      # Number of cells per simulation

        self.tmax = tmax                        # Total Time
        self.tau = 0                            # Tau Time
        self.time_stamp = time_stamp            # Time Stamp to save data
        self.division_time = division_time

        self._reference_division_time = (1/self.mu) * np.log(((self.birth_size + np.random.normal(loc=1, scale=0.05))/self.birth_size))
        self.reference_division_time = (1/self.mu) * np.log(((self.birth_size + np.random.normal(loc=1, scale=0.05))/self.birth_size))

        self.time_to_divide = 0

        # ___ Object Settings ___ #

        self.reaction_model = reaction_model    # Object of reaction model
        self._reaction_model_copy = copy.deepcopy(reaction_model)

        self.model_name = f'{model_name}.csv'

        # ___ Others ___ #
        self.burst = burst
        self.d_burst = d_burst
        self.tmp_chunk = []
        self.chunk_counter = 0
        self.chunk_lenght = chunk_lenght


    def simulate_division(self, model='adder'):
        
        if model == 'sizer':
            self.sizer()

        elif model == 'timer':
            pass

        elif model == 'adder':
            self.adder() # bs + np.random.normal(loc=1, scale=0.1)


    def simulate(self):
        """ 
            Performs Cesar Nieto et-al Gillespie algorithm
            for chemical reactions simulations.
        """
        schema = self.create_my_file(mode='classic')

        with open(self.model_name, mode='a', newline='') as f:

            writer = csv.DictWriter(f, fieldnames= schema, delimiter=',')
            
            for cell in tqdm(range(1, self.cells+1)):
             
                while self.time < self.tmax:
                    
                    possible_reactions = None

                    propensities = self.reaction_model.organize_propensities()
                    propenisties_list = list(propensities.values())

                    warnings.simplefilter("ignore")

                    calculated_tau = list(map(self.calculate_tau, propenisties_list))
                    reactions_dict = dict(enumerate(list(self.reaction_model.reactions.keys())))

                    self.tau = np.min(calculated_tau)
                    q = np.argmin(calculated_tau)       

                    if self.time + self.tau < self.reference_time:

                        self.reaction_model.solve_eqs(tau=self.tau)

                        possible_reactions = self.reaction_model.q[reactions_dict[q]]
                        self.react(reactions=possible_reactions)

                        self.time += self.tau
                        
                    elif self.time + self.tau > self.reference_time:
               
                        self.reaction_model.solve_eqs(tau=self.tau)
                        self.time = self.reference_time  
                        model = self.model_data(cell=cell)

                        if self.time_stamp != False:

                            if self.time >= self.time_stamp: 
                                writer.writerow(model)

                        else:
                            writer.writerow(model)

                        self.reference_time += self.sampling_time

                        model = None
                
                self.reaction_model = copy.deepcopy(self._reaction_model_copy)
                self.time = 0
                self.reference_time = 0
                
        f.close()


    def get_time_distribution(self, mode='single_switch', initial_state='off', species_name=None, threshold=0):

        """ 
            Performs Cesar Nieto et-al Gillespie algorithm
            for chemical reactions simulations.
        """

        schema = self.create_my_file(mode='time_distribution')

        with open(self.model_name, mode='a', newline='') as f:
            writer = csv.DictWriter(f, fieldnames= schema, delimiter=',')

            for cell in tqdm(range(self.cells)):
                
                on_counter = 0
                off_counter = 0
                for_dswitch = initial_state
                flag = True

                while flag:

                    possible_reactions = None

                    propensities = self.reaction_model.organize_propensities()
                    propenisties_list = list(propensities.values())

                    warnings.simplefilter("ignore")

                    calculated_tau = list(map(self.calculate_tau, propenisties_list))
                    reactions_dict = dict(enumerate(list(self.reaction_model.reactions.keys())))

                    self.tau = np.min(calculated_tau)
                    q = np.argmin(calculated_tau)       

                    self.reaction_model.solve_eqs(tau=self.tau)
                    possible_reactions = self.reaction_model.q[reactions_dict[q]]
                    self.react(reactions=possible_reactions)

                    self.time += self.tau

                    model = self.model_data_for_tau_time(cell=cell)

                    if mode == 'single_switch':
                        if initial_state == 'off':

                            if self.reaction_model.get_species(name=species_name) > threshold and on_counter == 0:
 
                                self.data_chunk(model=model, writer=writer)

                                on_counter = 1
                                flag = False

                        if initial_state == 'on':

                            if self.reaction_model.get_species(name=species_name) <= threshold and off_counter == 0:
                        
                                self.data_chunk(model=model, writer=writer)

                                off_counter = 1
                                flag = False
                            
                    elif mode == 'double_switch':

                        if for_dswitch == 'off':

                            if self.reaction_model.get_species(name=species_name) > threshold and on_counter == 0:
                    
                                self.data_chunk(model=model, writer=writer)

                                on_counter = 1
                                for_dswitch = 'on'

                        if for_dswitch == 'on':

                            if self.reaction_model.get_species(name=species_name) < threshold and off_counter == 0:
                        
                                self.data_chunk(model=model, writer=writer)

                                off_counter = 1
                                for_dswitch = 'off'
                            
                        if on_counter == 1 and off_counter == 1:
                            flag = False

        f.close()


    def sizer(self):

        """ 
            Performs Cesar Nieto et-al Gillespie algorithm
            for chemical reaction and division simulations.
        """

        schema = self.create_my_file(mode='division')

        with open(self.model_name, mode='a', newline='') as f:

            writer = csv.DictWriter(f, fieldnames= schema, delimiter='|')
            self.size = self.birth_size

            for cell in tqdm(range(1, self.cells+1)):
                
                while self.time < self.tmax:
                    reaction_q_copy = copy.deepcopy(self.reaction_model.q)
                    possible_reactions = None
                        
                    # _____Calculate Propensities_____ #
                    # - Discriminate between creation and degradataion

                    calculated_tau = self.calculate_sorted_tau()
                    calculated_tau['reference_division_time'] = self.reference_division_time

                    warnings.simplefilter("ignore")

                    calculated_tau_propensities = list(calculated_tau.values()) # Propensities
                    reactions_dict = dict(enumerate(list(calculated_tau.keys())))


                    self.tau = np.min(calculated_tau_propensities)
                    q = np.argmin(calculated_tau_propensities)
    
                    if self.time + self.tau < self.reference_time:
                        if reactions_dict[q] == 'reference_division_time':
                            reaction_q_copy['reference_division_time'] = {'reference_division_time': {'':''}}
                        
                        possible_reactions = reaction_q_copy[reactions_dict[q]]

                        self.size *= np.exp(self.mu*self.tau)
                        self.react(reactions=possible_reactions)     
                        self.time += self.tau               
                        
                    elif self.time + self.tau > self.reference_time:

                        self.reference_division_time -= (self.reference_time - self.time)
                        self.size *= np.exp(self.mu * (self.reference_time - self.time))
                        self.time = self.reference_time  

                        # ___Save Data___ # 
                        model = self.model_data_for_division(cell=cell, size=self.size, division_time=self.reference_division_time)

                        if self.time_stamp != False:

                            if self.time >= self.time_stamp: 
                                writer.writerow(model)

                        else:
                            writer.writerow(model)

                        self.reference_time += self.sampling_time

                self.reaction_model = copy.deepcopy(self._reaction_model_copy)
                self.time = 0
                self.reference_time = 0
                self.reference_division_time = 0


    def adder(self):

        """ 
            Performs Cesar Nieto et-al Gillespie algorithm
            for chemical reaction and division simulations.
        """

        schema = self.create_my_file(mode='division')

        with open(self.model_name, mode='a', newline='') as f:

            writer = csv.DictWriter(f, fieldnames= schema, delimiter='|')
            self.size = self.birth_size

            for cell in tqdm(range(1, self.cells+1)):
                while self.time < self.tmax:
                    reaction_q_copy = copy.deepcopy(self.reaction_model.q)
                    possible_reactions = None
                        
                    # _____Calculate Propensities_____ #
                    # - Discriminate between creation and degradataion

                    calculated_tau = self.calculate_sorted_tau()
                    calculated_tau['reference_division_time'] = self.reference_division_time

                    warnings.simplefilter("ignore")

                    calculated_tau_propensities = list(calculated_tau.values()) # Propensities
                    reactions_dict = dict(enumerate(list(calculated_tau.keys())))

                    self.tau = np.min(calculated_tau_propensities)
                    q = np.argmin(calculated_tau_propensities)
    
                    if self.time + self.tau < self.reference_time:
                        
                        self.reaction_model.solve_eqs(tau=self.tau)

                        if reactions_dict[q] == 'reference_division_time':
                            reaction_q_copy['reference_division_time'] = {'reference_division_time': {'':''}}
                        
                        possible_reactions = reaction_q_copy[reactions_dict[q]]

                        self.size *= np.exp(self.mu*self.tau)
                        self.react(reactions=possible_reactions)     
                        self.time += self.tau               
                        
                    elif self.time + self.tau > self.reference_time:

                        self.reaction_model.solve_eqs(tau=self.tau)

                        self.reference_division_time -= (self.reference_time - self.time)
                        self.size *= np.exp(self.mu * (self.reference_time - self.time))
                        self.time = self.reference_time  

                        # ___Save Data___ # 
                        model = self.model_data_for_division(cell=cell, size=self.size, division_time=self.reference_division_time)

                        if self.time_stamp != False:

                            if self.time >= self.time_stamp: 
                                writer.writerow(model)

                        else:
                            writer.writerow(model)

                        self.reference_time += self.sampling_time

                self.reaction_model = copy.deepcopy(self._reaction_model_copy)
                self.time = 0

                self.reference_time = 0
                self.reference_division_time = self._reference_division_time

                self.size = self.birth_size

        f.close()


    def react(self, reactions):

        if 'create' not in reactions.keys():
            pass

        elif 'create' in reactions.keys():
            create_species = reactions['create']
            [self.reaction_model.create(name=species) for species in create_species if species != None]


        if 'create_rna' not in reactions.keys():
            pass

        elif 'create_rna' in reactions.keys():
            create_mrna = reactions['create_rna']
            [self.reaction_model.create_rna(name=species) for species in create_mrna if species != None]


        if 'destroy' not in reactions.keys():
            pass

        elif 'destroy' in reactions.keys():
            destroy_species = reactions['destroy']
            [self.reaction_model.destroy(name=species) for species in destroy_species if species != None]


        if 'activate' not in reactions.keys():
            pass

        elif 'activate' in reactions.keys():
            activate_species = reactions['activate']
            [self.reaction_model.activate(name=species) for species in activate_species if species != None]


        if 'deactivate' not in reactions.keys():
            pass

        elif 'deactivate' in reactions.keys():
            deactivate_species = reactions['deactivate']
            [self.reaction_model.deactivate(name=species) for species in deactivate_species if species != None]


        if 'burst' not in reactions.keys():
            pass

        elif 'burst' in reactions.keys():
            burst_species = reactions['burst']
            [self.reaction_model.burst(name=species, fixed_burst=self.burst) for species in burst_species if species != None]


        if 'burst_degradation' not in reactions.keys():
            pass

        elif 'burst_degradation' in reactions.keys():
            burst_deg_species = reactions['burst_degradation']
            [self.reaction_model.burst_degradation(name=species, fix_d_burst=self.d_burst) for species in burst_deg_species if species != None]


        if 'reference_division_time' not in reactions.keys():
            self.reference_division_time -= self.tau
        
        elif 'reference_division_time' in reactions.keys():

            self.size /= 2
            self.reference_division_time = (1/self.mu) * np.log((self.size + np.random.normal(loc=1, scale=0.05))/self.size)   
            self.reaction_model.dilute_species() 
    

    def calculate_tau(self, propensity):
        warnings.simplefilter("ignore")
        p = -(1/propensity) * np.log(np.random.rand())
        return p


    def calculate_division_tau(self, propensity):
        warnings.simplefilter("ignore")
        p = (1/self.mu) * np.log(1 - (self.mu * np.log(np.random.rand())/(propensity * self.size)))
        return p


    def extract(self, key, iter_dict, values_to_get='keys'):
        
        if values_to_get == 'keys':
            return [propensity for propensity, value in iter_dict[key].items()]
        
        elif values_to_get == 'values':
            return [value for propensity, value in iter_dict[key].items()]
        
        else:
            raise KeyError


    def calculate_sorted_tau(self):

        propensity_type = self.reaction_model.propensity_type()

        reaction_type_I = ('create_rna', 'burst')
        reaction_type_II = ('destroy', 'activate', 'deactivate', 'create', 'burst_degradation')

        ultimate_propensities = {}

        create_propensity_type_I = None
        create_propensity_type_II = None

        for reaction_type in reaction_type_I:
                for key in propensity_type.keys():
                    if key == reaction_type:
                        propensity_values_type_I = self.extract(key=key, iter_dict=propensity_type, values_to_get='values')
                        propensity_to_create = list(map(self.calculate_division_tau, propensity_values_type_I))

                        propensity_keys_type_I = self.extract(key=key, iter_dict=propensity_type, values_to_get='keys')

                        create_propensity_type_I = dict(zip(propensity_keys_type_I, propensity_to_create))

                    else:
                        pass

        if create_propensity_type_I != None:
            ultimate_propensities.update(create_propensity_type_I)

        else:
            pass        


        for reaction_type in reaction_type_II:
            for key in propensity_type.keys():
            
                if key == reaction_type:
                    propensity_values_type_II = self.extract(key=key, iter_dict=propensity_type, values_to_get='values')
                    propensity_to_create = list(map(self.calculate_tau, propensity_values_type_II))

                    propensity_keys_type_II = self.extract(key=key, iter_dict=propensity_type, values_to_get='keys')

                    create_propensity_type_II = dict(zip(propensity_keys_type_II, propensity_to_create))

        warnings.simplefilter("ignore")
        
        if create_propensity_type_II != None:
            ultimate_propensities.update(create_propensity_type_II)
            
        else:
            pass    

        ultimate_propensities.update({'reference_division_time': 0})

        return ultimate_propensities


    def set_reaction(self, propensities, dart):
        sum_a = 0
        q = 0
        for i in range(len(propensities)):
            sum_a += propensities[i]
            if sum_a > dart:
                q = i
                break
        return q


    def create_my_file(self, mode='classic'):

        if mode == 'classic':
            schema = self.create_schema()

            try: 
                os.remove(self.model_name)
            except:
                pass
            try:
                with open(self.model_name, mode='x') as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=schema, delimiter=',')
                    writer.writeheader()
                    csvfile.close()
            except:
                print(f'File "{self.model_name}" already exist')
                pass

            return schema

        elif mode == 'division':
            
            schema = self.create_schema_for_division()

            try: 
                os.remove(self.model_name)
            except:
                pass
            try:
                with open(self.model_name, mode='x') as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=schema, delimiter=',')
                    writer.writeheader()
                    csvfile.close()
            except:
                print(f'File "{self.model_name}" already exist')
                pass

            return schema

        elif mode == 'time_distribution':
            
            schema = self.create_schema_for_tau_time()

            try: 
                os.remove(self.model_name)
            except:
                pass
            try:
                with open(self.model_name, mode='x') as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=schema, delimiter=',')
                    writer.writeheader()
                    csvfile.close()
            except:
                print(f'File "{self.model_name}" already exist')
                pass

            return schema


    def create_schema(self):
        schema = copy.deepcopy(self.reaction_model.species)
        schema = list(schema.keys())
        schema.append('time')
        schema.append('cell')
        return schema


    def model_data(self, cell):
        model = copy.deepcopy(self.reaction_model.species)
        model.update({'time': round(self.time, 2), 'cell': cell})
        return model


    def create_schema_for_division(self):
        schema = copy.deepcopy(self.reaction_model.species)
        schema = list(schema.keys())
        schema.append('time')
        schema.append('cell')
        schema.append('cell_size')
        schema.append('division_time')
        return schema


    def model_data_for_division(self, cell, size, division_time):
        model = copy.deepcopy(self.reaction_model.species)
        model.update({'time': round(self.time, 4), 'cell': cell, 'cell_size': round(size, 4), 'division_time': round(division_time, 4)})
        return model


    def model_data_for_tau_time(self, cell):
        model = copy.deepcopy(self.reaction_model.species)
        model.update({'time': self.time, 'cell': cell})
        return model


    def create_schema_for_tau_time(self):
        schema = copy.deepcopy(self.reaction_model.species)
        schema = list(schema.keys())
        schema.append('time')
        schema.append('cell')
        return schema


    def load_data(self):
        df = pd.read_csv(self.model_name, delimiter='|')
        return df


    def data_chunk(self, model, writer):
        self.tmp_chunk.append(model)
        self.chunk_counter += 1
        if self.chunk_counter == self.chunk_lenght:
            [writer.writerow(self.tmp_chunk[data]) for data in range(len(self.tmp_chunk))]
            self.chunk_counter = 0
            self.tmp_chunk = []

class ReactionModel:

    
    def __init__(self, reactions=dict,species=dict, propensities=dict, q=dict, parameters=dict, math_model=dict, exclude_from_division=[None]) -> None:
        self.reactions = reactions
        self.species = species 
        self.propensities = propensities
        self.q = q
        self.interaction_map = parameters
        self.math_model = math_model
        self.exclude_from_division = exclude_from_division
    

    def get_species(self, name=None):
        return self.species[name]


    def organize_propensities(self):
        
        complete_propensity = list(self.reactions.values())
        reaction_names = self.reactions.keys()

        species_copy = copy.copy(self.species)
        species_copy.update(self.propensities)

        species_calculus = [[species_copy[k] for k in complete_propensity[i]] for i in range(len(complete_propensity))]
        
        propensities_result = list(map(np.prod, species_calculus))
        return dict(zip(reaction_names, propensities_result))


    def sort_reactions(self):

        sorted_reactions = {key: self.reaction_type()[key] for key in self.reaction_type().keys() if len(self.reaction_type()[key]) >= 1}
        return sorted_reactions
        
    
    def propensity_type(self):

        propensities = self.organize_propensities()
        sorted_reactions = self.sort_reactions()

        for reaction_key, reaction_value in sorted_reactions.items():
            for reaction in reaction_value:
                sorted_reactions[reaction_key] = {key:propensities[key] for key in reaction_value}

        return sorted_reactions


    def reaction_type(self):

        storage = {
                    'create': [],
                    'create_rna': [],
                    'destroy':[],
                    'activate': [],
                    'deactivate': [], 
                    'burst': [], 
                    'burst_degradation': []
                    }

        for key, value in self.q.items():

            for reaction in storage.keys():

                storage[reaction].append(key) if reaction in list(value.keys()) else storage[reaction].append(None)
                {storage[k].remove(None) for k in storage.keys() if None in storage[k]}


        return storage


    def create(self, name=str):

        if self.species[name] >= 0:
            self.species[name] += 1

        elif self.species[name] < 0 :
            self.species[name] = 0
    

    def destroy(self, name=str):

        if self.species[name] > 0:
            self.species[name] -= 1

        elif self.species[name] <= 0 :
            self.species[name] = 0
    

    def activate(self, name=str):

        if self.species[name] >= 1 or self.species[name] <= 1:
            self.species[name] = 1
            

    def deactivate(self, name=str):

        if self.species[name] >= 0 or self.species[name] <= 0:
            self.species[name] = 0


    def burst(self, name=str, fixed_burst=15):

        self.species[name] += fixed_burst


    def burst_degradation(self, name=str, fix_d_burst=1.5):

        self.species[name] -= fix_d_burst


    def solve_eqs(self, tau=0):

        if self.math_model != dict:

            keys_to_save = list(self.math_model.keys())
            tmp_species = copy.deepcopy(self.species)
            tmp_species.update({'Tau': round(tau, 4)})

            values_list = [[[tmp_species[value] for value in list]] for list in list(self.interaction_map.values())]
            eqs_list = list(self.math_model.values())

            solved_eqs = [list(starmap(eqs_list[i], values_list[i])) for i in range(len(list(self.math_model.keys())))]
            values_to_save = [i[0] for i in solved_eqs]

            solution_to_save = {keys_to_save[idx]:round(values_to_save[idx], 4) for idx in range(len(keys_to_save))}
            
            for key in solution_to_save:
                self.species[key] = solution_to_save[key]

        else:
            pass

    
    def create_rna(self, name=str):

        if self.species[name] >= 0:
            self.species[name] += 1

        elif self.species[name] < 0 :
            self.species[name] = 0
    

    def dilute_species(self):

        species_to_dilute = list(self.species.keys())

        for name in species_to_dilute:

            if name not in self.exclude_from_division:
                diluted_species = np.random.binomial(self.species[name], 0.5) if self.species[name] > 0 else 0
                self.species[name] = diluted_species

            elif name in self.exclude_from_division:
                pass

if __name__ == '__main__':

    central_dogma = ReactionModel(
                species = {
                    'DNA': 1,
                    'mRNA': 0,
                    'Protein': 0
                },
                propensities = {
                    'transcription_propensity' : 0.3,
                    'translation_propensity': 10,
                    'mRNA_degradation' : np.log(2)/2,
                    'Protein_degradation' : np.log(2)/60
                },
                reactions = {
                    'transcription' : ['transcription_propensity', 'DNA'],
                    'translation' :   ['translation_propensity', 'mRNA'],
                    'mRNA_deg' :      ['mRNA_degradation', 'mRNA'],
                    'Protein_deg' :   ['Protein_degradation', 'Protein'],
                },
                q = {
                    'transcription' :   {'create_rna' : ['mRNA']},
                    'translation' :     {'create' : ['Protein']},
                    'mRNA_deg' :        {'destroy' : ['mRNA']},
                    'Protein_deg' :     {'destroy' : ['Protein']},
                }
)

    central_dogma_sim = Gillespie(reaction_model=central_dogma, model_name='central_dogma', tmax=50, sampling_time=1, cells=1)
    central_dogma_sim.simulate_division(model='adder')