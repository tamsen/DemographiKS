import math
import curve_fitting
from scipy.ndimage import gaussian_filter

class Ks_modeling_predictions:

    def __init__(self, config, bins):

        #bin_midpoints = [0.5 * (bins[i] + bins[i + 1]) for i in range(0, len(bins) - 1)]

        self.theoretical_ks_mean_now = config.mean_Ks_from_Tc + config.t_div_as_ks
        # theoretical_sigma_from_kingman_in_time= (2.0*config_used.ancestral_Ne)**-1 #1/K
        # theoretical_sigma_from_kingman_in_Ks= theoretical_sigma_from_kingman_in_time*config_used.Ks_per_YR
        # for an exponential, λ = 1/μ. And  σ =1/λ =μ. SO  σ =μ
        self.theoretical_Kingman_sigma_now=config.mean_Ks_from_Tc
        # theoretical exponential prediction
        K = config.mean_Ks_from_Tc ** -1
        bin_size = bins[1] - bins[0]
        popt = [config.num_genes * bin_size, config.t_div_as_ks, config.Ks_per_YR, K]
        self.travelling_kingman_ys=[curve_fitting.wgd_travelling_exponential(x, *popt) for x in bins]

        # theoretical gaussian prediction
        num_recombination_events_per_nuc=config.recombination_rate*config.WGD_time_Ge
        num_recombination_events_per_gene=num_recombination_events_per_nuc*config.gene_length_in_bases
        n = num_recombination_events_per_gene
        self.theoretical_Gaussian_sigma_now=self.theoretical_Kingman_sigma_now / math.sqrt(n)
        popt=[config.num_genes*bin_size,self.theoretical_ks_mean_now,self.theoretical_Gaussian_sigma_now]
        self.travelling_gaussian_ys = [curve_fitting.wgd_normal(x, *popt) for x in bins]

class Ks_modeling_fits:

    def __init__(self, slim_ks_by_gene, hist_Ns, bins):

        num_slim_genes = len(slim_ks_by_gene)
        self.mean_ks_now_from_slim = sum(slim_ks_by_gene) / num_slim_genes
        self.variance_from_slim = (1.0 / float(num_slim_genes)) * \
                             sum([(x - self.mean_ks_now_from_slim) ** 2 for x in slim_ks_by_gene])  # sigma_squared

        bin_midpoints = [0.5 * (bins[i] + bins[i + 1]) for i in range(0, len(bins) - 1)]
        self.gaussian_fit = curve_fitting.fit_curve_to_xs_and_ys(
            bin_midpoints, hist_Ns, curve_fitting.wgd_normal)

        self.gaussian_modified_results = curve_fitting.fit_curve_to_xs_and_ys(
            bin_midpoints, hist_Ns, curve_fitting.gaussian_modified_exponential)


def smooth_data(ys,sigma):

    smoothed_ys = gaussian_filter(ys, sigma=sigma)
    return smoothed_ys
