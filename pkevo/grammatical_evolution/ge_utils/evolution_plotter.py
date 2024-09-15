# Based on Robert Haas' "panakeias_garden" package (/panakeias_garden/pre_release_demos/metaheuristics.py)

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from config.file_storage_singleton import file_storage
import config.custom_logger as custom_logger

logger = custom_logger.get_logger(__name__, "ge_logs.log")

class EvolutionPlotter:
    """
    A class for visualizing evolutionary optimization progress.

    The `EvolutionPlotter` class is designed to visualize the progress of GE algorithm
    by displaying various statistics over generations, such as fitness scores and codon usage.

    Attributes:
        METRICS (int): The number of metrics used for statistical analysis, including
            minimum, maximum, mean, standard deviation and percentiles.

    """
    METRICS = 9 # min, max, mean, std, 10% percentile, 25% percentile, median, 75% percentile, 90% percentile

    def __init__(self, optimization_type, generations):
        """
        Initialize an EvolutionPlotter object

        Args:
            optimization_type (str): The type of optimization, either"max" or "min" (depends on selected Fitness function)
            generations (int): Number of forseen generations in the simulation (needed for size of data arrays)
        """
        self.optimization_type = optimization_type
        self.fitness_over_generations = np.zeros((self.METRICS, generations))
        self.codon_usage_over_generations = np.zeros((self.METRICS, generations))
        self.phenotypes_over_generations = np.zeros((3, generations))
    

    def _plot_metrics(self, measure_name, data, color, x, store_as_file):
        """
        Plot the results using Matplotlib.

        Args:
            measure_name (str): The name of the measure.
            data (np.ndarray): Measure values.
            color (str): Plot color.
            x (list): X-axis values (Number of the generation).
            store_as_file (bool): Store the plot as a file.
        """
        # Setup basic plot details
        fig = plt.figure()
        fig.set_figwidth(8)
        fig.suptitle(measure_name + " over generations")
        plt.xlabel('Generation')
        plt.ylabel(measure_name)
        # Create the actual plot inside a 1x1 subplot
        ax = plt.subplot(111)

        # Depending on optimization type, the trendline for best score should either follow the 'min' or the 'max' values
        if self.optimization_type == "max":
            z_best = np.polyfit(x, data[1], 1)
        else:
            z_best = np.polyfit(x, data[0], 1)

        if measure_name == "Fitness" and self.optimization_type == "max":
            ax.set_title("(Higher values are better)", fontdict={'fontsize':10})
        elif measure_name == "Fitness" and self.optimization_type == "min":
            ax.set_title("(Lower values are better)", fontdict={'fontsize':10})

        # Plot all 9 measures
        ax.plot(x, data[0], label="Min", color=color, linewidth=1, marker=".") # Min
        ax.plot(x, data[1], label="Max", color=color, linewidth=1, marker="2") # Max

        # Plots a line representing the Mean with errorbars (for the sd)
        ax.errorbar(x, data[2], yerr=data[3],
                    color='black', alpha=0.75, capthick=1, elinewidth=0.75,
                    label="Mean ± standard deviation")
        
        ax.plot(x, data[4], color='black', linestyle="solid", linewidth=0.0) # 10% Percentile
        ax.plot(x, data[5], color='black', linestyle="solid", linewidth=0.0) # 25% Percentile
        ax.plot(x, data[6], label="Median", color=color, linestyle="dotted", linewidth=1) # 50% Percentile, aka Median
        ax.plot(x, data[7], color='black', linestyle="solid", linewidth=0.0) # 75% Percentile
        ax.plot(x, data[8], color='black', linestyle="solid", linewidth=0.0) # 90% Percentile

        ax.fill_between(x, data[0], data[4], label="100% of the points", facecolor=color, alpha=0.075)  # within are all points
        ax.fill_between(x, data[8], data[1], facecolor=color, alpha=0.075)                              # within are all points
        ax.fill_between(x, data[4], data[5], label="80% of all points", facecolor=color, alpha=0.2)  # within are 80% of points
        ax.fill_between(x, data[7], data[8], facecolor=color, alpha=0.15)                            # within are 80% of points
        ax.fill_between(x, data[5], data[7], label="50% of all points", facecolor=color, alpha=0.3)  # within are 50% of points

        # Add trendlines for Mean and Best (based on https://www.statology.org/matplotlib-trendline/)
        # polyfit() returns coefficients of polynomial in ascending order -> first one defines the slope of the line
        p_best = np.poly1d(z_best)
        ax.plot(x, p_best(x), marker="x", color='orangered', alpha=0.9, linestyle="dashed", label="Trendline for Best")

        z_mean = np.polyfit(x, data[2], 1)
        p_mean = np.poly1d(z_mean)
        ax.plot(x, p_mean(x), color='orangered', alpha=0.9, linestyle="dashed", label="Trendline for Mean")

        # Add a legend
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True, ncol=1)
        
        # Remove the top and right bounding box lines
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)

        # Store the plot as a file or show it to the user directly
        if store_as_file:
            file_storage.save_plot(measure_name.replace(" ", "_"))
        else:
            plt.show()


    def _plot_phenotype_composition(self, data, x, store_as_file):   
        """
        Plot the results using Matplotlib.

        Args:
            data (np.ndarray): Measure values.
            x (list): X-axis values (Number of the generation).
            store_as_file (bool): Store the plot as a file.
        """
        # Setup basic plot details
        fig = plt.figure()
        fig.suptitle('Composition of Phenotypes Across Generations')
        plt.xlabel('Generation')
        plt.ylabel('Counts')
        # Create the actual plot inside a 1x1 subplot
        ax = plt.subplot(111)

        # Plot all 9 measures
        ax.bar(x, data[0], color='green', label='Unique Phenotypes')
        ax.bar(x, data[1], color='orange', bottom=data[0], label='Duplicated Phenotypes')
        ax.bar(x, data[2], color='red', bottom=data[0]+data[1], label='Invalid Phenotypes')

        # Add a legend
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True, ncol=1)
        
        # Remove the top and right bounding box lines
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)

        # Store the plot as a file or show it to the user directly
        if store_as_file:
            file_storage.save_plot('Phenotype_composition')
        else:
            plt.show()

    
    def add_data(self, fitness_stats, used_codons_stats, phenotype_stats, idx):
        """
        This function allows to iteratively add new data to the data arrays `self.fitness_over_generations`
        and `self.codon_usage_over_generations`. 

        Args:
            fitness_stats (list): Statistics about the achieved fitness scores
            used_codons_stats (list): Statistics about the used codons
            idx (int): The index of the data to add.
        """
        # Add the new data to the data arrays for each metric
        self.fitness_over_generations[:, idx] = fitness_stats
        self.codon_usage_over_generations[:, idx] = used_codons_stats
        self.phenotypes_over_generations[:, idx] = phenotype_stats


    def plot_evolution(self, iterations, store_as_file):
        """
        Initiates the actual plot generation after first preparing the plot data.

        Args:
            iterations (int): The number of iterations (=generations) performed by the evolutionary loop.
            store_as_file (bool): Store the plot as a file.
        """
        if iterations > 1:  
            # Our plots to be generated
            # Note: the data arrays might need to be cut after {iterations} in case the evolutionary loop ended prematurely
            fitness_plot = {"name": "Fitness", "data": self.fitness_over_generations[:, :iterations], "color": "mediumblue"}
            codon_usage_plot = {"name": "Codon usage", "data": self.codon_usage_over_generations[:, :iterations], "color": "forestgreen"}
            standard_plots = [fitness_plot, codon_usage_plot]

            # On the X-axis we will have the individual generations in consecutive order
            x_axis = range(iterations)

            with sns.axes_style("white"):
                for plot in standard_plots:
                    self._plot_metrics(plot["name"], plot["data"], plot["color"], x_axis, store_as_file)
                self._plot_phenotype_composition(self.phenotypes_over_generations[:, :iterations], x_axis, store_as_file)
        else:
            logger.warning("Evolution plotter needs at least 2 generations to create plots.")