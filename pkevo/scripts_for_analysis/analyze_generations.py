# Note: This script needs to be executed from root_dir/pkevo/ like this: python -m scripts_for_analysis.analyze_generations -f ../simulation_results/sim_2023-11-15_10-08-59/

import os
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import json
import argparse

def main(file_path):
    current_dir = os.path.dirname(os.path.abspath(__file__)) # should be "/pkevo/config"
    app_root_dir = os.path.join(current_dir, "..") # should be "/" 
    print(app_root_dir)
    file_path = os.path.join(app_root_dir, file_path, "snapshots.json")
    print(file_path)
    with open(file_path, 'r') as file:
        data = json.load(file)

    # Extract phenotypes for each generation
    phenotypes_per_generation = []
    result = dict()
    for item in data:
        generations = item['generation']
        generation_phenotypes = [gen['phenotype'] for gen in generations]
        phenotype_count = dict(Counter(generation_phenotypes))
        print(list(phenotype_count.values()))
        print(len(phenotype_count))
        print("----------")
        for key, value in phenotype_count.items():
            if value == 13:
                print(key)
        phenotypes_per_generation.append(generation_phenotypes)

    # Count occurrences of phenotypes in each generation
    occurrences_per_generation = []
    for generation in phenotypes_per_generation:
        occurrence_count = Counter(generation)
        occurrences_per_generation.append(occurrence_count)

    # Count unique and duplicate occurrences for each generation
    unique_counts_per_generation = []
    non_unique_counts_per_generation = []
    for generation_counts in occurrences_per_generation:
        unique_count = sum(count == 1 for count in generation_counts.values())
        duplicate_count = len(generation_counts) - unique_count
        unique_counts_per_generation.append(unique_count)
        non_unique_counts_per_generation.append(duplicate_count)

    # Plotting the bar graph
    generations = range(1, len(phenotypes_per_generation) + 1)
    plt.bar(generations, unique_counts_per_generation, label='Unique Phenotypes', color='blue')
    plt.bar(generations, non_unique_counts_per_generation, 
            label='Non unique Phenotypes', color='orange')

    plt.xlabel('Generation')
    plt.ylabel('Count')
    plt.title('Distribution of Unique and Duplicate Phenotypes per Generation')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-f', '--sim_folder', required=True, type=str,
                    help='Provide simulation folder name.')
    
    args = parser.parse_args()
    main(args.sim_folder)