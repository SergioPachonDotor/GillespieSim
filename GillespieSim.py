import numpy as np
from numpy.core.fromnumeric import sort, swapaxes
import pandas as pd
import copy
import warnings
import csv
import os
from tqdm import tqdm
from itertools import starmap

import matplotlib.pyplot as plt

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
                    birth_size=np.random.normal(loc=2/2, scale=0.1*(2/2))
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

        self.reference_division_time = (1/self.mu) * np.log(((self.birth_size + np.random.normal(loc=1, scale=0.05))/self.birth_size))
        self.time_to_divide = 0

        # ___ Object Settings ___ #

        self.reaction_model = reaction_model    # Object of reaction model
        self._reaction_model_copy = copy.deepcopy(reaction_model)

        self.model_name = f'{model_name}.csv'


    def simulate_division(self, model='sizer'):
        
        if model == 'sizer':
            self.sizer()

        elif model == 'timer':
            pass

        elif model == 'adder':
            pass # bs + np.random.normal(loc=1, scale=0.1)


    def simulate(self):
        """ 
            Performs Cesar Nieto et-al Gillespie algorithm
            for chemical reactions simulations.
        """
        schema = self.create_my_file()

        with open(self.model_name, mode='a', newline='') as f:

            writer = csv.DictWriter(f, fieldnames= schema, delimiter='|')
            
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


    def sizer(self):

        """ 
            Performs Cesar Nieto et-al Gillespie algorithm
            for chemical reaction and division simulations.
        """

        time_array = []
        mrna_array = []
        size_array = []

        mRNA=1
        mRNA_c=100
        mRNA_d=10
        
        # Tiempo actual de la simulación
        # Tiempo de referencia para el próximo muestreo

        schema = self.create_schema_for_division()

        with open(self.model_name, mode='a', newline='') as f:

            writer = csv.DictWriter(f, fieldnames= schema, delimiter='|')
            size = self.birth_size 

            while self.time < self.tmax:
        
                if mRNA == 0:

                    mRNA_synth = (1/self.mu) * np.log(1 - (self.mu * np.log(np.random.rand())/(mRNA_c * size)))

                    self.tau = np.min([mRNA_synth, self.reference_division_time])
                    q = np.argmin([mRNA_synth, self.reference_division_time])
                
                elif mRNA > 0:
                    
                    # _____Calculate Propensities_____ #
                    # - Discriminate between creation and degradataion

                    mRNA_synth = (1/self.mu) * np.log(1 - (self.mu * np.log(np.random.rand())/(mRNA_c * size)))
                    mRNA_deg = -(1/(mRNA_d * mRNA)) * np.log(np.random.rand())

                    self.tau = np.min([mRNA_synth,mRNA_deg, self.reference_division_time])
                    q = np.argmin([mRNA_synth,mRNA_deg, self.reference_division_time])

                if self.time + self.tau < self.reference_time:

                    size *= np.exp(self.mu*self.tau)

                    if q == 0:
                        mRNA += 1

                    elif q == 1:
                        mRNA -= 1
                    
                    elif q == 2:
                        size /= 2
                        self.reference_division_time = (1/self.mu) * np.log((size + np.random.normal(loc=1, scale=0.05))/size)
                        
                        # ____
                        mRNA = np.random.binomial(mRNA, 0.5)

                    if q != 2:
                        
                        self.reference_division_time -= self.tau

                    self.time += self.tau
                    
                    
                elif self.time + self.tau > self.reference_time:

                    self.reference_division_time -= (self.reference_time - self.time)
                    size *= np.exp(self.mu * (self.reference_time - self.time))

                    self.time = self.reference_time  

                    # ___Save Data___ # 
                    time_array.append(self.time)
                    mrna_array.append(mRNA)
                    size_array.append(size)


                    self.reference_time += self.sampling_time

            fig, ax1 = plt.subplots()

            ax2 = ax1.twinx()
            ax1.plot(time_array,size_array, '--',color='black', linewidth=1.5, label='size', alpha=0.7)
            ax2.plot(time_array,mrna_array,'ro', alpha=0.5, linewidth=0.5, ms=0.5, label='mRNA')

            ax1.set_xlabel('Time(h)')
            ax1.set_ylabel(r'Size (${\mu}$m)')
            ax2.set_ylabel('mRNA')
            plt.title('Implemented Adder')
            ax2.grid()
            plt.legend(loc=5)


    def react(self, reactions):


        if 'create' not in reactions.keys():
            pass

        elif 'create' in reactions.keys():
            create_species = reactions['create']
            [self.reaction_model.create(name=species) for species in create_species if species != None]


        if 'create_rna' not in reactions.keys():
            pass

        elif 'create_rna' in reactions.keys():
            create_species = reactions['create_rna']
            [self.reaction_model.create_rna(name=species) for species in create_species if species != None]


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
            create_species = reactions['burst']
            [self.reaction_model.burst(name=species) for species in create_species if species != None]


    def calculate_tau(self, propensity):
        return  -(1/propensity) * np.log(np.random.rand())


    def calculate_division_tau(self, propensity):
        return (1/self.mu) * np.log(1 - (self.mu * np.log(np.random.rand())/(propensity * self.size)))
        

    def extract(self, key, iter_dict, values_to_get='keys'):
        if values_to_get == 'keys':
            return [propensity for propensity, value in iter_dict[key].items()]
        
        elif values_to_get == 'values':
            return [value for propensity, value in iter_dict[key].items()]
        
        else:
            raise KeyError


    def calculate_sorted_tau(self):

        propensity_type = self.reaction_model.propensity_type()

        reaction_type_I = ('create', 'burst')
        reaction_type_II = ('destroy', 'activate', 'deactivate')

        ultimate_propensities = {}

        for reaction_type in reaction_type_I:
            for key in propensity_type.keys():

                if key == reaction_type:
                    propensity_values_type_I = self.extract(key=key, iter_dict=propensity_type, values_to_get='values')
                    propensity_to_create = list(map(self.calculate_division_tau, propensity_values_type_I))

                    propensity_keys_type_I = self.extract(key=key, iter_dict=propensity_type, values_to_get='keys')

                    create_propensity_type_I = dict(zip(propensity_keys_type_I, propensity_to_create))

                else:
                    pass

        for reaction_type in reaction_type_II:
            for key in propensity_type.keys():
            
                if key == reaction_type:
                    propensity_values_type_II = self.extract(key=key, iter_dict=propensity_type, values_to_get='values')
                    propensity_to_create = list(map(self.calculate_tau, propensity_values_type_II))

                    propensity_keys_type_II = self.extract(key=key, iter_dict=propensity_type, values_to_get='keys')

                    create_propensity_type_II = dict(zip(propensity_keys_type_II, propensity_to_create))

        ultimate_propensities.update(create_propensity_type_I)
        ultimate_propensities.update(create_propensity_type_II)

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
                    writer = csv.DictWriter(csvfile, fieldnames=schema, delimiter='|')
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
                    writer = csv.DictWriter(csvfile, fieldnames=schema, delimiter='|')
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
        schema.append('size')
        schema.append('division_time')
        return schema


    def model_data_for_division(self, cell, size, division_time):
        model = copy.deepcopy(self.reaction_model.species)
        model.update({'time': round(self.time, 2), 'cell': cell, 'size': size, 'division_time': division_time})
        return model


    def load_data(self):
        df = pd.read_csv(self.model_name, delimiter='|')
        return df


class ReactionModel:

    
    def __init__(self, reactions=dict,species=dict, propensities=dict, q=dict, interaction_map=dict, math_model=dict) -> None:
        self.reactions = reactions
        self.species = species 
        self.propensities = propensities
        self.q = q
        self.interaction_map = interaction_map
        self.math_model = math_model
    

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

        storage = {'create': [], 'destroy':[],'activate': [],'deactivate': [], 'burst': []}

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


    def burst(self, name=str):

        self.species[name] += 100


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


    def dilute_species(self):
        species_to_dilute = list(self.species.keys())
        for name in species_to_dilute:
            diluted_species = np.random.binomial(self.species[name], 0.5)
            self.species[name] = diluted_species
