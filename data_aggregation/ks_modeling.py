import math
import curve_fitting
from scipy.ndimage import gaussian_filter

class Ks_modeling_result:

    def __init__(self, config, bins):
       
       result=predict_Ks_from_config(config, bins)
       self.initial_kingman_as_ks=result[0]
       self.ks_model_exponential=result[1]
       self.ks_model_smoothed_exponential=result[2]
       self.ks_model_as_gaussian=result[3]


def fit_Ks_results(hist_results, bins):
    bin_midpoints = [0.5 * (bins[i] + bins[i + 1]) for i in range(0, len(bins) - 1)]
    gaussian_results = curve_fitting.fit_curve_to_xs_and_ys(
        bin_midpoints, hist_results, curve_fitting.wgd_gaussian)

    gaussian_modified_results = curve_fitting.fit_curve_to_xs_and_ys(
        bin_midpoints, hist_results, curve_fitting.gaussian_modified_exponential)
    return [gaussian_results,gaussian_modified_results]

def predict_Ks_from_config(config_used, bins):

    #synonymous_mutations_rate_only= config_used.mutation_rate / 1.2
    synonymous_mutations_rate_only = config_used.Ks_per_YR
    expected_Ks_peak_shift = config_used.DIV_time_Ge * synonymous_mutations_rate_only
    t_div_as_ks = config_used.DIV_time_Ge * config_used.Ks_per_YR
    bin_size = bins[1]-bins[0]
    print("config_used.mutation_rate " + str(config_used.mutation_rate))
    bin_size_in_time = bin_size / config_used.mutation_rate
    ks_for_one_generation= 1 * synonymous_mutations_rate_only
    two_Ne = 2.0 * config_used.bottleneck_Ne
    #bin_mid_points=[0.5*(bins[i]+bins[i+1]) for i in range(0,len(bins)-1)]

    kingman = [min(config_used.num_genes,
                   (bin_size_in_time * config_used.num_genes / two_Ne) * math.e ** (
                               (-1 * ((i / synonymous_mutations_rate_only) -1)) / two_Ne))
               for i in bins]
    Ks_model_exponential = []
    for b in bins:
        if b < expected_Ks_peak_shift-bin_size:
            Ks_model_exponential.append(0)
        else:
            Ks_model_exponential = Ks_model_exponential + kingman
            break

    ks_model_smoothed_exponential = smooth_data(Ks_model_exponential,(1/bin_size)/1000)
    #gaussian prediction:
    ks_model_as_gaussian = make_gaussian_prediction(bin_size, bins, config_used, expected_Ks_peak_shift, t_div_as_ks,
                                                    two_Ne)
    return kingman, Ks_model_exponential,ks_model_smoothed_exponential, ks_model_as_gaussian



def make_gaussian_prediction(bin_size, bins, config_used, expected_Ks_peak_shift, t_div_as_ks, two_Ne):
    total_ks_shift = config_used.mean_Ks_from_Tc + t_div_as_ks
    fit_fxn = curve_fitting.wgd_gaussian
    amp = config_used.num_genes * bin_size
    # The variance of an exponential distribution is 1/λ²
    # sigma= two_Ne*synonymous_mutations_rate_only / ( bin_size * 50)  # in Ks space
    # sigma = 0.002 #for Tdiv, Ne=1000, but Tdiv
    # seems to work for low Tdiv. Tdiv = 1000, Ne for 100,1000,5000
    sigma = two_Ne / 100000  # some fxn of RC (or, tending towards sigma of exponetial!)
    # mu=total_ks_shift
    mu = expected_Ks_peak_shift + config_used.mean_Ks_from_Tc
    popt = [amp, mu, sigma]
    ks_model_as_gaussian = [fit_fxn(b, *popt) for b in bins]
    return ks_model_as_gaussian


def smooth_data(ys,sigma):

    smoothed_ys = gaussian_filter(ys, sigma=sigma)
    return smoothed_ys
