import copy

from config.file_storage_singleton import file_storage

from .individual_encoder import IndividualEncoder

class GenerationSnapshots:
    """
    A class for storing and managing snapshots of generations created during the GE simulation.

    Attributes:
        generations (list): A list of generations, where each generation is a list of individuals.
    """

    def __init__(self):
        """
        Initialize an empty GenerationSnapshots object.
        """
        self.snapshots = []


    def add_generation(self, iteration, generation):
        """
        Add a new generation to the snapshots along with its iteration number.
        Args:
            iteration (int): The iteration number of the generation.
            generation: A list of individuals representing a single generation.
        """
        self.snapshots.append({
            'iteration': iteration,
            'generation': copy.copy(generation)
        })


    def get_generation(self, iteration):
        """
        Retrieve a specific generation by its iteration number.
        Args:
            iteration (int): The iteration number of the desired generation.
        Returns:
            tuple: A tuple containing the iteration number and the requested generation as a list of individuals.
        Raises:
            ValueError: If the specified iteration is not found.
        """
        for snapshot in self.snapshots:
            if snapshot['iteration'] == iteration:
                return snapshot['iteration'], snapshot['generation']
        raise ValueError("Iteration not found in snapshots")


    def get_last_generation(self):
        """
        Get the last (most recent) generation along with its iteration number.
        Returns:
            tuple: A tuple containing the iteration number and the last generation as a list of individuals.
        Raises:
            ValueError: If no generations are available.
        """
        if self.snapshots:
            iteration, generation = self.snapshots[-1]
            return iteration, generation
        else:
            raise ValueError("No generations available")


    def get_first_generation(self):
        """
        Get the first generation along with its iteration number.
        Returns:
            tuple: A tuple containing the iteration number and the first generation as a list of individuals.
        Raises:
            ValueError: If no generations are available.
        """
        if self.snapshots:
            iteration, generation = self.snapshots[0]
            return iteration, generation
        else:
            raise ValueError("No generations available")


    def get_generation_count(self):
        """
        Get the total number of stored generations.
        Returns:
            int: The number of stored generations.
        """
        return len(self.snapshots)
    

    def export_snapshots(self):
        """
        Export the stored snapshots of generations as a JSON file.

        This function serializes the stored snapshots of generations into a JSON file.
        Each snapshot includes the iteration number and the list of individuals in that generation.

        The exported JSON file is saved with the name 'snapshots.json'.
        """
        # Transform data, so that it can be handled as JSON
        snapshots_data = [{'iteration': snapshot['iteration'], 'generation': snapshot['generation']} for snapshot in self.snapshots]
        # Make the call to the File handler, which creates the JSON file
        file_storage.save_json_file('snapshots.json', snapshots_data, IndividualEncoder)

        self.visualize_pop_distribution()

    
    def visualize_pop_distribution(self):
        """
        """
        import matplotlib.pyplot as plt
        import numpy as np
        from collections import Counter
        
        # Initialize lists to store module counts for each snapshot
        all_modules_initial = []
        all_modules_final = []

        for snapshot in self.snapshots:
            unique_phenotypes = dict()
            init_pop_size = len(snapshot['generation'])

            modules = []

            for ind in snapshot['generation']:
                if ind.phenotype:
                    ind.build_pks_architecture()

                    modules.append(ind.architecture.module_count)

                    if ind.phenotype in unique_phenotypes:
                        unique_phenotypes[ind.phenotype] += 1
                    else:
                        unique_phenotypes[ind.phenotype] = 1
                else:
                    init_pop_size -= 1

            # Append module counts for the current snapshot
            if snapshot['iteration'] == 0:
                all_modules_initial.extend(modules)
            else:
                all_modules_final.extend(modules)

        # Count the occurrences of modules for each snapshot
        # counter_mod_initial = Counter(all_modules_initial)
        # counter_mod_final = Counter(all_modules_final)

        # Extract values for the x and y axes
        # x_values_mod_initial, y_values_mod_initial = zip(*sorted(counter_mod_initial.items()))
        # x_values_mod_final, y_values_mod_final = zip(*sorted(counter_mod_final.items()))
        
        # Create grouped bar plots
        plt.figure(figsize=(6, 5))
        values = [all_modules_initial, all_modules_final]
        bins = np.arange(1, max(max(all_modules_initial), max(all_modules_final)) + 2) - 0.5
        colors = ['#568ee9', '#e68a00']
        labels = ['Initial Population', 'Final Population']

        plt.hist(values, bins=bins, density=True, color=colors, label=labels, edgecolor='black')
        #plt.bar(np.array(x_values_mod_initial) - width, y_values_mod_initial, width=width, color='#568ee9', alpha=0.75, label='Initial Population')
        #plt.bar(np.array(x_values_mod_final) + width, y_values_mod_final, width=width, color='#e68a00', alpha=0.75, label='Final Population')

        plt.xlabel('Number of Modules')
        plt.ylabel('Frequency of PKSs')
        plt.title('PKS Distribution by Number of Modules')
        #plt.xticks(np.arange(min(min(x_values_mod_initial), min(x_values_mod_final)), max(max(x_values_mod_initial), max(x_values_mod_final))+1, 1))
        plt.legend()

        plt.tight_layout()

        file_storage.save_plot("population_distribution")
        #plt.show()