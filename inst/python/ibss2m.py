import numpy as np
import jax.numpy as jnp
from inst.python.ser import fit_tilted_ser, fit_uvb_ser
from inst.python.ser import initialize_re, initialize_ser_params

def pluck(dicts: dict, k) -> np.ndarray:
    """Pull key from dictionary of dictionaries

    Args:
        dicts (dict): a dictionary of dictionaries, each entry has the same keys
        k (str): a key in the dictionary

    Returns:
        np.ndarray : array stacking the results of all dictionaries
    """
    if isinstance(k, list) and len(k) > 1:
        k0 = k[0]
        k = k[1:]
        peeled = {key: value[k0] for key, value in dicts.items()}
        res = pluck(peeled, k)
    else:
        if isinstance(k, list):
            k = k[0]
        res = np.array([val[k] for _, val in dicts.items()])
    return res


def stack(dicts: dict):
    """Stack multiple dictionaries with shared set of keys.

    Args:
        dicts (dict): dictionary of dictionaries

    Returns:
        dict: a dictionary with the same keys as each element of dicts
    """

    # TODO: can we do this recursively?
    # good if we want to organize the output of SER function
    example_key = list(dicts.keys())[0]
    keys = dicts[example_key].keys()
    return {k : pluck(dicts, k) for k in keys}


def compute_pips(alpha: jnp.array):
    pips = 1 - np.exp(np.nansum(np.log(1 - alpha), axis=0))
    return pips


def add_re(re, ser, keep_2m=True):
    new_re = dict(
        mu = re['mu'] + ser['summary']['psi'],
        mu2 = (re['mu2'] + ser['summary']['psi2'] + 2 * re['mu'] * ser['summary']['psi'])
    )
    if not keep_2m:
        new_re['mu2'] = new_re['mu']**2
    return new_re


def subtract_re(re, ser, keep_2m=True):
    new_re = dict(
        mu = re['mu'] - ser['summary']['psi'],
        mu2 = (re['mu2'] - ser['summary']['psi2'] - 2 * (re['mu'] - ser['summary']['psi']) * ser['summary']['psi'])
    )
    if not keep_2m:
        new_re['mu2'] = new_re['mu']**2
    return new_re


def check_offset(re: dict):
    var = re['mu2'] - re['mu']**2
    return(jnp.all(var >= 0))


def make_ibss2m(ser_fun):
    
    def ibss2m_iter(data, sers, re, control):
        """
        complete 1 loop over SERs in IBSS2M
        """
        L = len(sers)
        n, p = data['X'].shape
        keep_2m = control.get('keep_2m', True)

        diff = 0
        new_sers = dict()
        for l in range(L):
            # remove the effect of old SER
            ser_l = sers[l]
            re = subtract_re(re, ser_l, keep_2m)

            # fit new SER
            params = ser_l['params']
            
            # params = initialize_ser_params(n, p, ser_l['tau0'])
            new_ser_l = ser_fun(data, re, params, {})

            # add back effect, record SER
            re = add_re(re, new_ser_l, keep_2m)
            new_sers[l] = new_ser_l
            diff = diff + jnp.sum((ser_l['summary']['psi'] - new_ser_l['summary']['psi'])**2)
        return new_sers, re, diff

    def ibss2m_driver(data: dict , params: dict, control: dict):
        """Fit IBSS Second Moments. Inputs are expected to be complete,
        handling default settings deferred to a wrapper function.

        Args:
            data (dict): dictionary with data, keys X and y
            params (dict): L, prior_variance, prior_weights
            control (dict): maxit estimate_intercept, estimate_prior_variance

        Returns:
            dict: a dictionary with fit variational parameters and posterior summaries
        """
        n, p = data['X'].shape

        L = params['L']
        pi = params['pi']
        tau0 = params['tau0']

        # ignore second moment-- reverts to IBSS
        keep_2m = control.get('keep_2m', True)

        # first iteration, initialize SERs and add to re
        # tau0 = 1 / params['prior_variance'] 
        re = initialize_re(n)
        sers = dict()
        for l in range(L):
            ser_params = initialize_ser_params(n, p, tau0)
            ser_params['pi'] = pi  # use prior weights specified
            ser_l = ser_fun(data, re, ser_params, {})# fit SER
            re = add_re(re, ser_l, keep_2m)  # add SER to random effect
            sers[l] = ser_l  # store fit SER

        # iterate ibss2m until convergence
        diff_history = []
        elbo_history = []
        kl_history = []
        lbf_history = []
        for i in range(control['maxit']):
            sers, re, diff = ibss2m_iter(data, sers, re, control)
            diff_history.append(diff)
            elbo_history.append(pluck(sers, ['track', 'elbo_ser']))
            kl_history.append(pluck(sers, ['track', 'kl_ser']))
            lbf_history.append(pluck(sers, ['summary', 'lbf']))
            if (diff < control['tol']): 
                break 

        # 0. stack up SER results
        # elbo_hsitory and kl_history give the elbo and KL of each SER
        # which can be used to get the ELBO for the full model
        track = dict(
            diff = diff,
            converged = diff < control['tol'],
            it=i,
            diff_history=np.array(diff_history),
            elbo_history=np.array(elbo_history),
            kl_history=np.array(kl_history),
            lbf_history=np.array(lbf_history)
        )
        res = dict(sers=sers, re=re, track=track)
        return res

    def ibss2m_jax(X, y, L = 10,
        prior_variance = 1.,
        prior_weights = None,
        tol = 1e-3,
        maxit = 100,
        estimate_intercept = True,
        estimate_prior_variance = False,
        keep_2m = True):

        data = dict(X=np.array(X), y=np.array(y))

        if prior_weights is None:
            p = X.shape[1]
            prior_weights = jnp.ones(p)/p

        params = dict(
            L=L,
            tau0 = 1/prior_variance,
            pi = prior_weights
        )
        
        control = dict(
            tol=tol,
            maxit=maxit,
            estimate_intercept=estimate_intercept,
            estimate_prior_variance=estimate_prior_variance,
            keep_2m = keep_2m
        )

        fit = ibss2m_driver(data, params, control)

        alpha = pluck(fit['sers'], ['params', 'alpha'])
        mu = pluck(fit['sers'], ['params', 'mu'])
        tau = pluck(fit['sers'], ['params', 'tau'])
        lbf = pluck(fit['sers'], ['summary', 'lbf'])
        pip = compute_pips(alpha)
        psi = dict(
            mu = np.array(fit['re']['mu']),
            mu2 = np.array(fit['re']['mu2'])
        )
        [fit['sers'][l]['params'].pop('xi') for l in fit['sers']];

        susie = dict(
            sers = fit['sers'],
            alpha = alpha,
            mu = mu,
            var = 1/tau,
            lbf = lbf,
            pip = pip,
            psi = psi,
            track = fit['track']
        )
        return susie

    return ibss2m_jax
    
ibss2m_tilted = make_ibss2m(fit_tilted_ser)
ibss2m_uvb = make_ibss2m(fit_uvb_ser)
