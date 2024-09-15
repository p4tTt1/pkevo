# Credits: Based on the LigandScout wrapper class from Robert Haas' unified_cheminformatics package.

import subprocess
import os
import re
from numbers import Number
from rdkit.Chem import PandasTools 

from config.pkevo_config import GeneralConfig
from config.file_storage_singleton import file_storage
import config.custom_logger as custom_logger

logger = custom_logger.get_logger(__name__, log_file='ligand_logs.log')

class LigandScoutUtils():
    """
    Utility class to access necessary LigandScout (or short = `lig`) CLI tools
    for performing fitness evaluation of generated polyketides against
    target receptors.

    Access is accomplished via Python subprocesses.

    This class is implemented as a Singleton, ensuring that there is only one
    instance of the class, allowing consistent access to the LigandScout CLI
    tools.

    Methods:
        get_ligand_scout_version(): Determine the version of the installed
            toolkit and return "unknown" in case of failure.
        get_pharmacophore_score_list: Retrieve a list of pharmacophore
            scores for a given list of SMILES strings using the specified
            pharmacophore and optional parameters.

    """
    _instance = None

    def __new__(cls): 
        """
        Create a new instance of the class as a Singleton.

        Returns:
            LigandScoutUtils: The Singleton instance of the class.
        """
        lig_path = GeneralConfig.LIGANDSCOUT_PATH
        max_score = 0.0

        # ensure path contains a trailing path separator
        if not lig_path.endswith(os.path.sep):
            lig_path += os.path.sep

        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance.lig_path = lig_path
            cls._instance.max_score = max_score
        return cls._instance


    def _command_builder(self, program_name, arguments=None, keyword_arguments=None, assignment_symbol='='):
        """
        Build a command string for executing a command-line program with arguments and keyword arguments.

        Args:
            program_name (str): The name of the command-line program to be executed.
            arguments (list, optional): List of arguments to be included in the command string.
            keyword_arguments (dict, optional): Dictionary of keyword arguments and their values
                to be included in the command string.
            assignment_symbol (str, optional): The symbol used to assign values to keyword arguments
                in the command string. Default is '='.

        Returns:
            str: The complete command string that can be executed in a subprocess.
        """
        # Argument processing
        if arguments is None:
            arguments = []
        if keyword_arguments is None:
            keyword_arguments = []

        # Transformation
        cmd = program_name
        for arg in arguments:
            cmd += ' ' + str(arg)
        for key, value in keyword_arguments.items():
            if value is not None:
                cmd += ' ' + str(key) + assignment_symbol + str(value)

        logger.info(f"Executing command: {cmd}")
        return cmd

    
    def _run_lig_cli_command(self, command, timeout=None):
        """
        Execute a command in a subprocess using the LigandScout command-line interface (CLI).

        This function runs the specified command in a subprocess, capturing its standard output
        and standard error streams, and returning the standard output if the command is successful.

        Args:
            command (str): The command to be executed, as a single string.
            timeout (float, optional): Maximum time (in seconds) for the command to run before
                raising a TimeoutError. Default is None (no timeout).

        Returns:
            str: The standard output of the executed command if it succeeds, or None if there
            is an error or if the command times out.
        """
        try:
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True, timeout=timeout)
            if result.returncode == 0:
                return result.stdout.strip()
            else:
                logger.warning(f"There seems to be a problem with LigandScout: {result.stderr.strip()}")
                return None
        except subprocess.TimeoutExpired as timeout_err:
            logger.error("TimeoutError: There seems to be a problem with the connection the LigandScout. Please check that you can access LigandScout from your computer.")
            return None
        except Exception as e:
            logger.error("Error:", str(e))
            return None


    def _check_scores_validity(self, name_list, score_list):
        """
        Check if a list of scores corresponds to a list of names and contains valid numbers.

        This function verifies whether the number of names in the name list matches the number of scores
        in the score list and whether all the scores in the score list are valid numbers (instances of Number).

        Args:
            name_list (list): A list of names to be associated with the scores.
            score_list (list): A list of scores to be validated.

        Raises:
            ValueError: If the number of names and scores are not equal, or if any score is not a valid number.
        """
        if len(name_list) != len(score_list):
            raise ValueError('Number of names ({}) and number of scores ({}) are not '
                             'equal.'.format(len(name_list), len(score_list)))
        for score in score_list:
            if not isinstance(score, Number):
                raise ValueError('Score is not a valid number: {}'.format(score))


    ## aka run_idbgen()
    def _smi_file_to_conformation_database(
            self, smi_filepath, ldb_filepath, log_file=None, set_memory=None, confgen_type=None,
            multi_conf=None, num_confs=None, name_property=None, omega_licence=None,
            hydrophob_threshold=None, error=None, update_keep_old=None, update_drop_old=None,
            allow_duplicates=None, overwrite_existing=None, num_processes=None, slave_memory=None,
            max_molweight=None, exec_mode=None, num_cores=None, first_index=None, last_index=None,
            index_file=None, index_db=None, start_positions=None, host_list_file=None,
            comma_sep_hosts=None, kill_workers=None, network_interface=None, verbosity=None):
        """
        Execute the LigandScout 'idbgen' command to generate a conformation database from a SMILES file.

        Args:
            smi_filepath (str): Path to the input SMILES file.
            ldb_filepath (str): Path to the output LigandScout database (LDB) file.
            Note: Remaining arguments are LigandScout-specific. See its manual or help page for 
                `idbgen` tool.

        Returns:
            str: Path to the generated LigandScout database (LDB) file.
        """
        # Transformation
        command = self._command_builder(
            program_name = self.lig_path + 'idbgen',
            keyword_arguments = {
                '--input': smi_filepath,
                '--output': ldb_filepath,
                '--log-file': log_file,
                '--set-memory': set_memory,
                '--confgen-type': confgen_type,
                '--multi-conf': multi_conf,
                '--num-confs': num_confs,
                '--name-property': name_property,
                '--omega-license': omega_licence,
                '--hydrophob_threshold': hydrophob_threshold,
                '--error': error,
                '--update-keep-old': update_keep_old,
                '--update-drop-old': update_drop_old,
                '--allow-duplicates': allow_duplicates,
                '--overwrite-existing': overwrite_existing,
                '--num-processes': num_processes,
                '--slave-memory': slave_memory,
                '--max-molweight': max_molweight,
                '--exec-mode': exec_mode,
                '--num-cores': num_cores,
                '--first-index': first_index,
                '--last-index': last_index,
                '--index-file': index_file,
                '--index-db': index_db,
                '--start-positions': start_positions,
                '--host-list-file': host_list_file,
                '--comma-sep-hosts': comma_sep_hosts,
                '--kill-workers': kill_workers,
                '--network-interface': network_interface,
                '--verbosity': verbosity
            },
            assignment_symbol = ' ')

        result = self._run_lig_cli_command(command)
        #print(result)
        return ldb_filepath


    ## aka _run_iscreen()
    def _conformation_database_to_pharmacophore_sdf_file(
            self, ldb_filepath, pharmacophore_filepath, sdf_filepath,
            log_file=None, simple_progress=None, set_memory=None, conformation_match_mode=None,
            disable_exclusion_volumes=None, fragment_screening_mode=None, allow_omit=None,
            min_required=None, regex_file=None, first_index=None, last_index=None, num_cores=None,
            scoring_function=None, boolean_expression=None, roc_curve=None, timeout_cl=None,
            exec_mode=None, host_list_file=None, comma_sep_hosts=None, kill_workers=None,
            network_interface=None, verbosity=None):
        """
        Execute the LigandScout 'iscreen' command to perform pharmacophore screening on a conformation database.

        Args:
            ldb_filepath (str): Path to the LigandScout database (LDB) file.
            pharmacophore_filepath (str): Path to the pharmacophore file used for screening.
            sdf_filepath (str): Path to the output SDF file to store the screening results.
            scoring_function (str, optional): Scoring function used for screening.
            Note: Remaining arguments are LigandScout-specific. See its manual or help page for 
                `idbgen` tool.

        Returns:
            str: Path to the output SDF file containing pharmacophore screening results.

        """
        # Transformation
        command = self._command_builder(
            program_name= self.lig_path + 'iscreen',
            keyword_arguments = {
                '--query': pharmacophore_filepath,
                '--database': ldb_filepath,
                '--output': sdf_filepath,
                '--log-file': log_file,
                '--simple-progress': simple_progress,
                '--set-memory': set_memory,
                '--conformation-match-mode': conformation_match_mode,
                '--disable-exclusion-volumes': disable_exclusion_volumes,
                '--fragment-screening-mode': fragment_screening_mode,
                '--allow-omit': allow_omit,
                '--min-required': min_required,
                '--regex-file': regex_file,
                '--first-index': first_index,
                '--last-index': last_index,
                '--num-cores': num_cores,
                '--scoring-function': scoring_function,
                '--boolean-expression': boolean_expression,
                '--roc-curve': roc_curve,
                '--timeout-cl': timeout_cl,
                '--exec-mode': exec_mode,
                '--host-list-file': host_list_file,
                '--comma-sep-hosts': comma_sep_hosts,
                '--kill-workers': kill_workers,
                '--network-interface': network_interface,
                '--verbosity': verbosity
            },
            assignment_symbol = ' ')

        result = self._run_lig_cli_command(command)
        #print(result)
        return sdf_filepath


    def _smi_list_to_smi_file(self, smi_list, smi_file_path, name_list):
        """
        Write a list of SMILES strings to a .smi file with associated names.

        Args:
            smi_list (list of str): A list of SMILES strings representing chemical structures.
            smi_file_path (str): The path to the output .smi file.
            name_list (list of str): A list of names corresponding to the SMILES strings.

        """
        with open(smi_file_path, 'w') as file_handle:
            for smi, name in zip(smi_list, name_list):
                file_handle.write(smi+'\t'+name+'\n')


    def _process_sdf_result_file(self, filepath):
        """
        Read a pharmacophore screening SDF file and extract names and associated scores.

        Args:
            filepath (str): Path to the input SDF file generated from pharmacophore screening.

        Returns:
            Tuple[list, list]: A tuple containing two lists:
                - A list of names associated with the scored molecules.
                - A list of scores extracted from the SDF file.

        Raises:
            ValueError: If the SDF file structure does not match the expected format.
        """
        names = []
        scores = []

        # Use RDKit to load and read the SDF file
        sdf_df = PandasTools.LoadSDF(filepath)
        if not sdf_df.empty:
            names = sdf_df['ID'].tolist() # IDs are the 'names' (= row index) we gave the molecules when creating the .smi file
            scores = sdf_df['Score'].astype(float).tolist() # Score contains the values achieved by the selected Pharmacophore Fit function

            # Check validity of extracted data
            self._check_scores_validity(names, scores)

            # Check what the highest score from this run was and update overall max score if necessary 
            max_score_from_run = max(scores)
            if max_score_from_run > self.max_score:
                self.max_score = max_score_from_run
                # Additionally extract the best scored conformation and store it in the simulation results
                filtered_row = sdf_df[sdf_df['Score'] == str(max_score_from_run)] # Note: in the original df all values are strings
                #file_storage.save_dataframe_as_sdf(filtered_row)
                file_storage.move_sdf_results(filepath)

        return names, scores


    def get_ligand_scout_version(self):
        """
        Determine the version of the installed LigandScout toolkit.

        Returns:
            str: The LigandScout toolkit version, or "unknown" in case of failure.
        """
        tool_name = "iscreen"
        lig_command = " --help"
        full_command = self.lig_path + tool_name + lig_command
        
        try:
            result = self._run_lig_cli_command(full_command, timeout=5)
            pattern = 'LigandScout V(.+)'
            match = re.search(pattern, result)
            version = match.group(1)
        except Exception:
            version = 'unknown'
        return version
    

    def get_pharmacophore_score_list(
            self, smiles_list, pharmacophore_file, min_config=True, scoring_function='relative', 
            log_file=None, set_memory=None, confgen_type=None, multi_conf=None, num_confs=None,
            name_property=None, omega_licence=None, hydrophob_threshold=None, error=None,
            update_keep_old=None, update_drop_old=None, allow_duplicates=None,
            overwrite_existing=None, num_processes=None, slave_memory=None, max_molweight=None,
            exec_mode=None, num_cores=None, first_index=None, last_index=None, index_file=None,
            index_db=None, start_positions=None, host_list_file=None, comma_sep_hosts=None,
            kill_workers=None, network_interface=None, verbosity=None,
            simple_progress=None, conformation_match_mode=None, disable_exclusion_volumes=None,
            fragment_screening_mode=None, allow_omit=None, min_required=None, regex_file=None,
            boolean_expression=None, roc_curve=None, timeout_cl=None):
        """
        Retrieve a list of pharmacophore scores for a given list of SMILES strings.

        Args:
            smiles_list (list of str): A list of SMILES strings representing
                chemical structures to be evaluated.
            pharmacophore_file (str): The filename of the pharmacophore file to be used
                for screening.
            min_config (bool, optional): If True, use minimal configuration for idbgen
                and iscreen. If False, provide additional configuration options.
            scoring_function (str, optional): The scoring function to be used for
                pharmacophore screening.
            Note: Remaining arguments are LigandScout-specific. See its manual or help page for 
                details.

        Returns:
            list of float: A list of pharmacophore scores corresponding to the input SMILES
            structures. The order of scores in the list matches the order of input SMILES.
        """

        # Note to self: first we need to run idbgen and after that we can execute iscreen!
        try: 
            # Get the current directory where this script is located
            current_dir = os.path.dirname(os.path.abspath(__file__))
            pharmacophore_file_path = os.path.join(current_dir, '..', 'chemistry', 'chemistry_model', 'pharmacophores' ,pharmacophore_file)

            # create temporary folder and .smi file with all smiles for an individual
            tmp_dir_path = file_storage.create_tmp_dir()
            smi_file_path = os.path.join(tmp_dir_path, 'molecules.smi')
            ldb_file_path = os.path.join(tmp_dir_path, 'conformers.ldb')
            sdf_file_path = os.path.join(tmp_dir_path, 'conformers_scored.sdf')
    
            # for each molecule use its row index as name
            name_list_being_positions = [str(i) for i in range(len(smiles_list))]
            # create the temporary .smi file
            self._smi_list_to_smi_file(smiles_list, smi_file_path, name_list_being_positions)

            # Generate conformers (by executing `idbgen` with either minimal required arguments or all possible)
            if min_config:
                self._smi_file_to_conformation_database(smi_filepath=smi_file_path, ldb_filepath=ldb_file_path)
            else:
                self._smi_file_to_conformation_database(
                    smi_filepath=smi_file_path, ldb_filepath=ldb_file_path,
                    log_file=log_file, set_memory=set_memory, confgen_type=confgen_type,
                    multi_conf=multi_conf, num_confs=num_confs, name_property=name_property,
                    omega_licence=omega_licence, hydrophob_threshold=hydrophob_threshold, error=error,
                    update_keep_old=update_keep_old, update_drop_old=update_drop_old,
                    allow_duplicates=allow_duplicates, overwrite_existing=overwrite_existing,
                    num_processes=num_processes, slave_memory=slave_memory,
                    max_molweight=max_molweight,
                    exec_mode=exec_mode, num_cores=num_cores, first_index=first_index,
                    last_index=last_index, index_file=index_file, index_db=index_db,
                    start_positions=start_positions, host_list_file=host_list_file,
                    comma_sep_hosts=comma_sep_hosts, kill_workers=kill_workers,
                    network_interface=network_interface, verbosity=verbosity)

            # Pharmacophore screening (by executing `iscreen` with either minimal required arguments or all possible)
            if min_config:
                self._conformation_database_to_pharmacophore_sdf_file(
                    ldb_filepath=ldb_file_path, pharmacophore_filepath=pharmacophore_file_path,
                    sdf_filepath=sdf_file_path, scoring_function=scoring_function)
            else:
                self._conformation_database_to_pharmacophore_sdf_file(
                    ldb_filepath=ldb_file_path, pharmacophore_filepath=pharmacophore_file_path,
                    sdf_filepath=sdf_file_path, scoring_function=scoring_function,
                    log_file=log_file, simple_progress=simple_progress, set_memory=set_memory,
                    conformation_match_mode=conformation_match_mode,
                    disable_exclusion_volumes=disable_exclusion_volumes,
                    fragment_screening_mode=fragment_screening_mode, allow_omit=allow_omit,
                    min_required=min_required, regex_file=regex_file, first_index=first_index,
                    last_index=last_index, num_cores=num_cores,
                    boolean_expression=boolean_expression, roc_curve=roc_curve, timeout_cl=timeout_cl,
                    exec_mode=exec_mode, host_list_file=host_list_file,
                    comma_sep_hosts=comma_sep_hosts, kill_workers=kill_workers,
                    network_interface=network_interface, verbosity=verbosity)

            # Initiliaze Scoring results list with zeros
            scores_all = [0.0 for _ in range(len(name_list_being_positions))]
            # Get results from the SDF file
            names_being_positions, scores_calculated = self._process_sdf_result_file(sdf_file_path)
            # Map results from SDF into scoring results list
            for name, score in zip(names_being_positions, scores_calculated):
                idx = int(name)
                scores_all[idx] = score
            return scores_all

        except Exception as err:
            raise err
        finally:
            file_storage.delete_tmp_dir()


# Singleton instance of LigandScoutUtils (Note: import instance, not class!)
ligand_scout_utils = LigandScoutUtils()
