import math
from scipy.optimize import curve_fit
from scipy.stats import norm,exponnorm

#https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.exponnorm.html
def gaussian_modified_exponential(x, Amp, K, loc, scale):
    return Amp * exponnorm.pdf(x, K, loc, scale)

def coal(x, N):
    k = 1 / N
    return k * math.exp(k*(1 - x))

def linear(x, m, b):
    return m * x + b
def wgd_gaussian(x, amp, mu, sig):
    return amp * norm.pdf(x, mu, sig)

#    kingman = [min(num_slim_genes,
#                   (bin_size_in_time * num_slim_genes / two_Ne) * math.e ** (
#                               (-1 * (i / config_used.mutation_rate)) / two_Ne))
#               for i in bins]

def travelling_exp(x,Amp,K,O):
    if x <= O:
        return 0

    return Amp*math.exp(-K*(x-O))

def my_decay(x, Ao,k, r):
    return Ao*((1-r)**(k*x))


def fit_curve_to_xs_and_ys(xs_for_wgd, ys_for_wgd, fit_fxn,p0=False):

    try:
        if p0:
            popt, pcov = curve_fit(fit_fxn, xs_for_wgd, ys_for_wgd, p0=p0)
        else:
            popt, pcov = curve_fit(fit_fxn, xs_for_wgd, ys_for_wgd)

    except Exception as inst:
        print(type(inst))  # the exception type
        print(inst.args)  # arguments stored in .args
        print(inst)  # __str__ allows args to be printed directly,
        return False, xs_for_wgd, False

    fit_curve_ys = [fit_fxn(x, *popt) for x in xs_for_wgd]
    RMSE_to_sum = [(fit_curve_ys[i] - ys_for_wgd[i]) * (fit_curve_ys[i] - ys_for_wgd[i]) for i in
                   range(0, len(xs_for_wgd))]
    RMSE = math.sqrt(sum(RMSE_to_sum) / len(RMSE_to_sum))

    # rms2 = mean_squared_error(ys_for_wgd, fit_curve_ys, squared=False)
    # chi2 = do_chi2(fit_curve_ys, ys_for_wgd)


    return fit_curve_ys, xs_for_wgd, popt