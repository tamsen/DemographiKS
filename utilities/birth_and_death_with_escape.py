import math
import random

import numpy
import numpy as np
from matplotlib import pyplot as plt

def gene_birth_death_with_escape_pdf(
        step_size,
        decay_constant,
        escape_rate,
        rseed,
        include_debugging_plots):

    #this is the range to normalize the pdf
    max_ks = 4.0
    xs = np.arange(0, max_ks, step_size)

    exp_decay_ys=[math.exp(-1*decay_constant*x) for x in xs]
    ys_that_have_died=[1.0-math.exp(-1*decay_constant*x) for x in xs]
    ys_that_escape=[y*escape_rate for y in ys_that_have_died]
    ys_remaining=[exp_decay_ys[idx]+ys_that_escape[idx] for idx in range(0,len(exp_decay_ys))]

    normalization_factor=1.0/(step_size*sum(ys_remaining))
    normalized_pdf=[normalization_factor*y for y in ys_remaining]
    
    if include_debugging_plots:
        plt.close()
        plt.plot(xs, exp_decay_ys, label="SSDs not yet decayed", color='k')
        plt.plot(xs, ys_that_have_died, label="SSDs should be dead", color='b')
        plt.plot(xs, ys_that_escape, label="SSDs that escape", color='c')
        plt.plot(xs, ys_remaining, label="SSDs_remaining_over_time", color='r')
        plt.legend()
        out_png = "/home/tamsen/Data/DemographiKS_output_from_mesx/birth_death_model.png"
        plt.savefig(out_png)
        plt.close()
        
    return normalized_pdf, xs


def draw_SSDs_from_pdf(normalized_pdf, xs, num_genes_needed, rseed, include_debugging_plots):

        population = xs
        probabilities = normalized_pdf
        #random.seed(rseed)
        random_draws = random.choices(
            population, weights=probabilities, k= num_genes_needed)

        if include_debugging_plots:
            hist_ys_real, bins, patches = plt.hist(random_draws, bins=40, facecolor='b', alpha=0.25,
                                                    density=True, label="observation")
            out_hist_png = "/home/tamsen/Data/DemographiKS_output_from_mesx/bd_exp_hist.png"
            plt.savefig(out_hist_png)
            plt.close()

        return random_draws