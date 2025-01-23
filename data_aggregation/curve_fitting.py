import math
from scipy.optimize import curve_fit
from scipy.stats import norm,exponnorm

#https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.exponnorm.html
def gaussian_modified_exponential(x, Amp, K, loc, scale):
    return Amp * exponnorm.pdf(x, K, loc, scale)


def wgd_normal(x, amp, mu, sig):
    scale=(2.0*math.pi*sig*sig)**(-0.5)
    exponent=(x-mu)*(x-mu)/(2.0*sig*sig)
    return  amp * scale * math.e**(-1*exponent)

def wgd_travelling_exponential(x, amp, loc_of_maximum, K):

    if x <= loc_of_maximum:
        return 0
    else:
        result= amp * K * math.e ** (-K * ((x - loc_of_maximum)))
        return result


def mix_of_travelling_exponential_with_gaussian(x, amp_exp, loc_of_maximum,
                                                K,  mu, sig):

    if x <= loc_of_maximum:
        return 0
    else:
        exp_result= amp_exp * K * math.e ** (-K * ((x - loc_of_maximum)))
        exponent = (x - mu) * (x - mu) / (2.0 * sig * sig)
        #scale=(2.0*math.pi*sig*sig)**(-0.5)
        scale=1.0
        gaus_result = scale * math.e ** (-1 * exponent)
        return exp_result*gaus_result


def travelling_kingman_old(x, two_Ne, Ks_per_YR, bin_size_in_time, num_genes, tdiv_in_ks):

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