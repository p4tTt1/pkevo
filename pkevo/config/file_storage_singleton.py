import json
import os
import sys
import datetime
import shutil
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools


class FileStorage:
    _instance = None

    def __new__(cls):
        current_dir = os.path.dirname(os.path.abspath(__file__)) # should be "/pkevo/config"
        app_root_dir = os.path.join(current_dir, "../..") # should be "/" 
        simulation_base_dir = os.path.join(app_root_dir, "simulation_results")
        logs_dir = os.path.join(app_root_dir, "logs")
        tmp_dir = os.path.join(app_root_dir, "tmp_folder")

        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance.simulation_base_dir = simulation_base_dir
            cls._instance.logs_dir = logs_dir
            cls._instance.tmp_dir = tmp_dir
            #cls._instance.create_simulation_folder()
            #cls._instance.create_logs_folder()
        return cls._instance


    def create_simulation_folder(self):
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        simulation_folder_name = f"sim_{timestamp}"
        # define the 'simulation_folder_path' attribute inside the function and not during initialization
        # that way we can skip file writing requests if there is no folder for the simulation
        self.simulation_folder_path = os.path.join(self.simulation_base_dir, simulation_folder_name)
        if not os.path.exists(self.simulation_folder_path):
            os.makedirs(self.simulation_folder_path)


    def create_logs_folder(self):
        if not os.path.exists(self.logs_dir):
            os.makedirs(self.logs_dir)


    def create_tmp_dir(self):
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
        return os.path.abspath(self.tmp_dir)


    def delete_tmp_dir(self):
        if os.path.exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)

    def save_dataframe_as_sdf(self, dataframe):
        if hasattr(self, 'simulation_folder_path'):
            new_file_path = os.path.join(self.simulation_folder_path, 'best_conformation.sdf')
            PandasTools.WriteSDF(df=dataframe, out=new_file_path, properties=list(dataframe.columns))
    def move_sdf_results(self, sdf_file_path):
        if hasattr(self, 'simulation_folder_path'):
            new_file_path = os.path.join(self.simulation_folder_path, 'conformers_scored.sdf')
            shutil.move(sdf_file_path, new_file_path)

    def save_text_file(self, file_name, content):
        if hasattr(self, 'simulation_folder_path'):
            file_path = os.path.join(self.simulation_folder_path, file_name)
            with open(file_path, "w" , encoding="utf-8") as file:
                file.write(content)

    
    def save_json_file(self, file_name, content, encoder=None):
        if hasattr(self, 'simulation_folder_path'):
            file_path = os.path.join(self.simulation_folder_path, file_name)
            with open(file_path, 'w') as json_file:
                json.dump(content, json_file, cls=encoder)
        
    
    def prepare_individuals(self):
        """
        Read data from a JSON file and prepare a list of individuals.

        Returns:
            list: A list of potential Individuals represented by genome and phenotype.
        """
        file_path = "config/starting_individuals.json"
        with open(file_path, 'r') as file:
            data = json.load(file)

        individuals = []
        for entry in data:
            genome = entry.get("Genome", [])
            phenotype = entry.get("Phenotype", "")
            individual = [genome, phenotype]
            individuals.append(individual)

        return individuals
    

    def save_parameters(self, evolution_config, search_config):
        if hasattr(self, 'simulation_folder_path'):
            parameter_file_path = os.path.join(self.simulation_folder_path, "meta_information.txt")
            command = " ".join(sys.argv)

            with open(parameter_file_path, "w", encoding="utf-8") as parameter_file:
                parameter_file.write(f"Executed command: {command}\n")
                parameter_file.write("=================================\n\n")
                evo_params = [(a, getattr(evolution_config, a)) for a in dir(evolution_config) if not a.startswith('__')]
                parameter_file.write("Values used for Evolution Config:\n")
                parameter_file.write("=================================\n")
                for name, value in evo_params:
                    parameter_file.write(f"{name}: {value}\n")
                parameter_file.write("\n")

                search_params = [(a, getattr(search_config, a)) for a in dir(search_config) if not a.startswith('__')]
                parameter_file.write("Values used for Search Config:\n")
                parameter_file.write("==============================\n")
                for name, value in search_params:
                    parameter_file.write(f"{name}: {value}\n")
                parameter_file.write("\n")


    def save_best_ind(self, ind):
        if hasattr(self, 'simulation_folder_path'):
            output_file_path = os.path.join(self.simulation_folder_path, "best_polyketide.txt")
            with open(output_file_path, "w", encoding="utf-8") as output_file:
                output_file.write("Best individual:\n")
                output_file.write(ind.__str__(verbosity='vvv') + "\n")


    def save_plot(self, plot_name):
        if hasattr(self, 'simulation_folder_path'):
            import matplotlib.pyplot as plt
            plt.savefig(os.path.join(self.simulation_folder_path, plot_name + ".png"), bbox_inches="tight")
            plt.close()

    
    def save_molecule(self, smiles_str):
        if hasattr(self, 'simulation_folder_path'):
            result_mol = Chem.MolFromSmiles(smiles_str)

            # check shouldn't be necessary, since the fitness function already checked, but better safe than sorry 
            if result_mol is not None:
                output_file_path = os.path.join(self.simulation_folder_path, "result_molecule.svg")
                Draw.MolToFile(result_mol, output_file_path, imageType="svg", size=(300, 300), kekulize=True, wedgeBonds=True)
            

# Singleton instance of FileStorage (Note: import instance, not class!)
file_storage = FileStorage()
