import numpy as np
import pandas as pd
import copy
import warnings
import csv
import os
from tqdm import tqdm
from itertools import starmap

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
                    division_size=2
                    ):

        # __ Parameters ___ #
        
        # Time Parameters
        self.sampling_time = sampling_time      # Sampling Time
        self.time = 0                           # Iitial Time
        self.reference_time = 0                 # Reference Time
        self.tmax = tmax                        # Total Time
        self.tau = 0                            # Tau Time
        self.time_stamp = time_stamp            # Time Stamp to save data

        # Size Parameters
        self.mu = np.log(2)/division_time       # Mu parameter for Sizer Model
        self.division_size = division_size
        self.birth_size = np.random.normal(loc=self.division_size/2, scale=0.1*(self.division_size/2))

        # Cells
        self.cells = cells                      # Number of cells per simulation

        # ___ Object Settings ___ #

        self.reaction_model = reaction_model    # Object of reaction model
        self._reaction_model_copy = copy.deepcopy(reaction_model)

        self.model_name = f'{model_name}.csv'


    def simulate(self):
        """ 
            Performs Cesar Nieto et-al Gillespie algorithm
            for chemical reactions simulations.
        """
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

        with open(self.model_name, mode='a', newline='') as f:

            writer = csv.DictWriter(f, fieldnames= schema, delimiter='|')
            
            for cell in tqdm(range(1, self.cells+1)):
                # progress_bar = tqdm(total=self.tmax)
                while self.time < self.tmax:

                    propensities = self.reaction_model.calculate_propensities()
                    propenisties_list = list(propensities.values())
                    possible_reactions = None

                    warnings.simplefilter("ignore")

                    # reaction_times = lambda propensity, random_num: -(1/propensity) * np.log(random_num)
                    random_array = [np.random.rand() for k in range(len(propenisties_list))]

                    calculated_tau = list(map(self.calculate_tau, propenisties_list, random_array))
                    reactions_dict = dict(enumerate(list(self.reaction_model.show_reactions().keys())))
                    
                    self.tau = np.min(calculated_tau)     # tau time
                    q = np.argmin(calculated_tau)         # Reaction that occurs

                    if self.time + self.tau < self.reference_time:

                        # Intergate Eqs
                        self.reaction_model.solve_eqs(tau=self.tau)
                        #
                        possible_reactions = self.reaction_model.show_q()[reactions_dict[q]]

                        if 'create' not in possible_reactions.keys():
                            pass

                        elif 'create' in possible_reactions.keys():
                            create_species = possible_reactions['create']
                            [self.reaction_model.create(name=species) for species in create_species if species != None]


                        if 'destroy' not in possible_reactions.keys():
                            pass

                        elif 'destroy' in possible_reactions.keys():
                            destroy_species = possible_reactions['destroy']
                            [self.reaction_model.destroy(name=species) for species in destroy_species if species != None]


                        if 'activate' not in possible_reactions.keys():
                            pass

                        elif 'activate' in possible_reactions.keys():
                            activate_species = possible_reactions['activate']
                            [self.reaction_model.activate(name=species) for species in activate_species if species != None]


                        if 'deactivate' not in possible_reactions.keys():
                            pass

                        elif 'deactivate' in possible_reactions.keys():
                            deactivate_species = possible_reactions['deactivate']
                            [self.reaction_model.deactivate(name=species) for species in deactivate_species if species != None]    

                        self.time += self.tau
                        
                    else:
                        #
                        self.reaction_model.solve_eqs(tau=self.tau)
                        #
                        self.time = self.reference_time  
                        model = self.model_data(cell=cell)
                        if self.time_stamp != False:
                            if self.time >= self.time_stamp: 
                                writer.writerow(model)
                        else:
                            writer.writerow(model)
                        # progress_bar.update(self.sampling_time)

                        self.reference_time += self.sampling_time
                        model = None
                
                self.reaction_model = copy.deepcopy(self._reaction_model_copy)
                self.time = 0
                self.reference_time = 0
                # progress_bar.close()
        f.close()


    def simulate_gillespie(self):
        """ 
            Performs Gillespie-Nieto simualtion algorithm
            for chemical reactions simulations.
        """
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

        with open(self.model_name, mode='a', newline='') as f:

            writer = csv.DictWriter(f, fieldnames= schema, delimiter='|')
            for cell in tqdm(range(1, self.cells+1)):
                # progress_bar = tqdm(total=self.tmax)
                while self.time < self.tmax:

                    propensities = self.reaction_model.calculate_propensities()

                    possible_reactions = None

                    warnings.simplefilter("ignore")
                    a_total = [np.sum(list(propensities.values()))]
                    random_number = [np.random.rand()]

                    self.tau = list(map(self.calculate_tau, a_total, random_number))[0]

                    reactions_dict = dict(enumerate(list(self.reaction_model.show_reactions().keys())))

                    dart = a_total[0] * np.random.rand()    # Reaction that occurs

                    
                    if self.time + self.tau < self.reference_time:
                        propensities_to_use = list(propensities.values())
                        q = self.set_reaction(propensities=propensities_to_use, dart=dart)

                        # Intergate Eqs
                        self.reaction_model.solve_eqs(tau=self.tau)
                        #
                        possible_reactions = self.reaction_model.show_q()[reactions_dict[q]]

                        if 'create' not in possible_reactions.keys():
                            pass

                        elif 'create' in possible_reactions.keys():
                            create_species = possible_reactions['create']
                            [self.reaction_model.create(name=species) for species in create_species if species != None]


                        if 'destroy' not in possible_reactions.keys():
                            pass

                        elif 'destroy' in possible_reactions.keys():
                            destroy_species = possible_reactions['destroy']
                            [self.reaction_model.destroy(name=species) for species in destroy_species if species != None]


                        if 'activate' not in possible_reactions.keys():
                            pass

                        elif 'activate' in possible_reactions.keys():
                            activate_species = possible_reactions['activate']
                            [self.reaction_model.activate(name=species) for species in activate_species if species != None]


                        if 'deactivate' not in possible_reactions.keys():
                            pass

                        elif 'deactivate' in possible_reactions.keys():
                            deactivate_species = possible_reactions['deactivate']
                            [self.reaction_model.deactivate(name=species) for species in deactivate_species if species != None]         

                        self.time += self.tau
                        
                    else:
                        #
                        self.reaction_model.solve_eqs(tau=self.tau)
                        # print(f'{self.time + self.tau} > {self.reference_time}')
                        #
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


    def simulate_division(self):
        """ 
            Performs Gillespie-Nieto simualtion algorithm
            for chemical reactions simulations.
        """
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

        with open(self.model_name, mode='a', newline='') as f:

            writer = csv.DictWriter(f, fieldnames= schema, delimiter='|')
            
            for cell in tqdm(range(1, self.cells+1)):
                # progress_bar = tqdm(total=self.tmax)
                while self.time < self.tmax:

                    propensities = self.reaction_model.calculate_propensities()

                    possible_reactions = None

                    warnings.simplefilter("ignore")
                    a_total = [np.sum(list(propensities.values()))]
                    random_number = [np.random.rand()]

                    self.tau = list(map(self.calculate_tau, a_total, random_number))[0]

                    reactions_dict = dict(enumerate(list(self.reaction_model.show_reactions().keys())))

                    dart = a_total[0] * np.random.rand()    # Reaction that occurs

                    division_time = (1/self.mu) * np.log(self.division_size/self.birth_size)

                    # reaction_first = 'update_data' if self.reference_time < division_time else 'save_data'
                    
                    if self.time + self.tau < self.reference_time:

                        propensities_to_use = list(propensities.values())
                        q = self.set_reaction(propensities=propensities_to_use, dart=dart)

                        # Intergate Eqs
                        self.reaction_model.solve_eqs(tau=self.tau)
                        #
                        possible_reactions = self.reaction_model.show_q()[reactions_dict[q]]

                        if 'create' not in possible_reactions.keys():
                            pass

                        elif 'create' in possible_reactions.keys():
                            create_species = possible_reactions['create']
                            [self.reaction_model.create(name=species) for species in create_species if species != None]


                        if 'destroy' not in possible_reactions.keys():
                            pass

                        elif 'destroy' in possible_reactions.keys():
                            destroy_species = possible_reactions['destroy']
                            [self.reaction_model.destroy(name=species) for species in destroy_species if species != None]


                        if 'activate' not in possible_reactions.keys():
                            pass

                        elif 'activate' in possible_reactions.keys():
                            activate_species = possible_reactions['activate']
                            [self.reaction_model.activate(name=species) for species in activate_species if species != None]


                        if 'deactivate' not in possible_reactions.keys():
                            pass

                        elif 'deactivate' in possible_reactions.keys():
                            deactivate_species = possible_reactions['deactivate']
                            [self.reaction_model.deactivate(name=species) for species in deactivate_species if species != None]         

                        self.time += self.tau
                        
                    else:
                        #
                        self.reaction_model.solve_eqs(tau=self.tau)
                        # print(f'{self.time + self.tau} > {self.reference_time}')
                        #
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


    def calculate_tau(self, popensity, random_num):
        return -(1/popensity) * np.log(random_num)


    def set_reaction(self, propensities, dart):
        sum_a = 0
        q = 0
        for i in range(len(propensities)):
            sum_a += propensities[i]
            if sum_a > dart:
                q = i
                break
        return q
                
        
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


    def create_schema_for_division(self):
        schema = copy.deepcopy(self.reaction_model.show_species())
        schema = list(schema.keys())
        schema.append('time')
        schema.append('cell')
        schema.append('size')
        return schema


    def model_data_for_division(self, cell, size):
        model = copy.deepcopy(self.reaction_model.show_species())
        model.update({'time': round(self.time, 2), 'cell': cell, 'size': size})
        return model


    def load_data(self):
        df = pd.read_csv(self.model_name, delimiter='|')
        return df

    
    def sizer(self):
        pass


class Bacteria:

    def __init__(self) -> None:
        pass


class ReactionModel:

    
    def __init__(self, reactions=dict,species=dict, propensities=dict, q=dict, interaction_map=dict, math_model=dict) -> None:
        self.reactions = reactions
        self.species = species 
        self.propensities = propensities
        self.q = q
        self.interaction_map = interaction_map
        self.math_model = math_model
    

    def show_species(self):
        return self.species


    def show_reactions(self):
        return self.reactions


    def show_propensities(self):
        return self.propensities
    

    def show_q(self):
        return self.q


    def show_eqs(self):
        return self.math_model


    def show_solve(self):
        return self.solve


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


    def divide(self, name=str):
        diluted_species = np.random.binomial(self.species[name], 0.5)
        self.species[name] = diluted_species

