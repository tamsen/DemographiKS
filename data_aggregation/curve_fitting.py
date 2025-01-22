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

def wgd_exponential(x, amp, loc_of_maximum, K):

    if x <= loc_of_maximum:
        return 0
    else:
        result= amp * K * math.e ** (-K * ((x - loc_of_maximum)))
        return result


def travelling_kingman(x, two_Ne, Ks_per_YR, bin_size_in_time, num_genes, tdiv_in_ks):

    if x <= tdiv_in_ks:
        return 0
    else:
        scalar=bin_size_in_time * num_genes / two_Ne
        result= scalar * math.e ** ((-1 * (((x-tdiv_in_ks) /Ks_per_YR) -1)) / two_Ne)
        return min(num_genes,result)


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