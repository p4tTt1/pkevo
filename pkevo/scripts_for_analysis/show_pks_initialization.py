# Note: This script needs to be executed from root_dir/pkevo/ like this: python -m scripts_for_analysis.show_pks_initialization

import os
import matplotlib.pyplot as plt
import numpy as np

from grammatical_evolution.grammatical_evolution import GrammaticalEvolution

def ge_instance():
    # Get the current directory where this script is located
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the path to the grammar file
    grammar_file_path = os.path.join(current_dir, '..', 'grammatical_evolution', 'grammars', 'pks_simplified.bnf')

    # Create the GrammaticalEvolution object
    ge = GrammaticalEvolution(grammar_file_path, 
                              fitness_function_name="StringDistance",
                              verbose=False, 
                              store_results=False, 
                              use_starting_inds=False)
    return ge


def main(ge_instance: GrammaticalEvolution):
    init_pop_size = 1000
    pop = ge_instance.initialize_population(size=init_pop_size, starting_inds=None)
    unique_phenotypes = dict()

    # Initialize lists to store the number of domains and number of modules
    domains = []
    modules = []
    counter_3_mods = 0
    counter_4_mods = 0

    # Trigger fitness evaluation, so that the counter for unique and total Individuals are activated
    ge_instance.perform_evaluation(pop)

    for ind in pop:
        if ind.phenotype:
            ind.build_pks_architecture() # with this we now get the amount of modules and domains for each individual

            if ind.architecture.module_count == 3:
                counter_3_mods +=1
            elif ind.architecture.module_count == 4:
                counter_4_mods +=1

            domains.append(ind.architecture.domain_count)
            modules.append(ind.architecture.module_count)

            if ind.phenotype in unique_phenotypes:
                unique_phenotypes[ind.phenotype] += 1
            else:
                unique_phenotypes[ind.phenotype] = 1
        else:
            init_pop_size -= 1

    # Count the number of unique strings
    num_unique_strings = len(unique_phenotypes)

    avg_modules = sum(modules) / init_pop_size
    print(f"{avg_modules=}")
    avg_domains = sum(domains) / init_pop_size
    print(f"{avg_domains=}")
    print(f"{counter_3_mods=}")
    print(f"{counter_4_mods=}")

    ## Plots data
    # Count the occurrences of each value in the list
    counts_dom = {domain: domains.count(domain) for domain in set(domains)}
    counts_mod = {module: modules.count(module) for module in set(modules)}

    # Extract values for the x and y axes
    x_values_dom = list(counts_dom.keys())
    y_values_dom = list(counts_dom.values())
    x_values_mod = list(counts_mod.keys())
    y_values_mod = list(counts_mod.values())

    # Create bar charts
    plt.figure(figsize=(12, 5))

    plt.subplot(1, 2, 1)
    plt.bar(x_values_dom, y_values_dom, color='blue', alpha=0.75)
    plt.xlabel('Number of Domains')
    plt.ylabel('Number of PKSs')
    plt.title('Bar Plot for Number of Domains')
    # Customize x-axis tick marks and labels
    plt.xticks(range(min(x_values_dom), max(x_values_dom)+1, 6))

    plt.subplot(1, 2, 2)
    plt.bar(x_values_mod, y_values_mod, color='green', alpha=0.75)
    plt.xlabel('Number of Modules')
    plt.ylabel('Number of PKSs')
    plt.title('Bar Plot for Number of Modules')
    plt.xticks(range(min(x_values_mod), max(x_values_mod)+1, 1))

    plt.tight_layout()

    # Show the histograms
    plt.show()

    # Find strings that occur more than once
    duplicate_strings = {k: v for k, v in unique_phenotypes.items() if v > 1}

    print(f"Number of generated phenotypes: {init_pop_size}")
    print(f"Number of unique phenotypes: {num_unique_strings}")
    #print("Duplicate strings and their occurrences:")
    #for string, count in duplicate_strings.items():
    #    print(f"'{string}': {count} times")



if __name__ == "__main__":
    ge = ge_instance()
    main(ge)