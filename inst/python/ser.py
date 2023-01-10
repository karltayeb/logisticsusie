import jax.numpy as jnp
import numpy as np
from utils import categorical_kl, map_nested_dicts
from tilted import tilted_univariate_lr
from univariate_vb import fit_univariate_vb

from jax import jit, vmap
from jax.scipy.special import logsumexp

def initialize_re(n: int, mu: float = 0, mu2: float = 0):
    """Initialize empty random effects dicitonary

    Args:
        n (int): number of obsrvations

    Returns:
        dict: dictionary with mu and mu2
    """
    return dict(
        mu = jnp.ones(n) * mu,
        mu2 = jnp.ones(n) * mu2
    )

def initialize_univariate_params(n: int, tau0: float):
    mu_init = 0
    tau_init = 1
    xi_init = np.ones(n) * 1e-3
    delta_init = 0
    return dict(
        mu = mu_init,
        tau = tau_init,
        xi = xi_init,
        delta = delta_init,
        tau0 = tau0
    )

def initialize_ser_params(n: int, p: int, tau0: float):
    mu_init = np.zeros(p)
    tau_init = np.ones(p)
    xi_init = np.ones((n, p)) * 1e-3
    delta_init = np.zeros(p)
    tau0 = np.ones(p) * tau0
    pi = np.ones(p) / p
    alpha = pi
    return dict(
        alpha = alpha,
        mu = mu_init,
        tau = tau_init,
        xi = xi_init,
        delta = delta_init,
        tau0 = tau0,
        pi = pi
    )

def make_ser(fit_fun):
    """Create SER from univariate logistic regression fit.

    Args:
        fit_fun (function): a function for fitting unviariate logistic regression.

    Returns:
        function: function for fitting SER
    """
    # vectorize the fit function
    fit_fun_vec = jit(vmap( 
        fit_fun, 
        in_axes=({'X': 1, 'y': None}, {'mu': None, 'mu2':None}, {'mu': 0, 'tau': 0, 'xi': 1, 'delta': 0, 'tau0': 0}),
        out_axes={'params': {'mu': 0, 'tau': 0, 'xi': 1, 'delta': 0, 'tau0': 0}, 'track': {'elbo': 0, 'kl': 0, 'diff': 0, 'it': 0}}
    ))

    # function for fitting null model
    def tilted_null_model(y: jnp.array, re: dict):
        """Fit null model where the effect is exactly 0

        Args:
            y (jnp.array): n vector of binary observations
            re (dict): random effects term a 'mu' and 'mu2'
        """
        n = y.size
        X = y[::-1]
        data = dict(X=X, y=y)
        params = initialize_univariate_params(n, 1e10)
        null_model = fit_fun(data, re, params)
        return null_model
    
    # make SER function
    def fit_ser_jax(data: dict, re: dict, params: dict, control : dict):
        """
        data: dict with X and y
        re: dict with mu and mu2, first and second moments of random effect per observation
        params: dict with fitting parameters
        control: dict estimate_prior_variance
        """
        # 0: compute UVB for each observation
        pi = params.pop('pi')  # the univariate VBs don't know about pi
        alpha = params.pop('alpha')
        res = fit_fun_vec(data, re, params)

        # 1. compute pips (alpha) by normalizing elbos
        logits = res['track']['elbo'] + jnp.log(pi)  # unnormalized
        alpha = jnp.exp(logits - logsumexp(logits))
        elbo = np.sum(res['track']['elbo'] * alpha) - categorical_kl(alpha, pi)
        kl = np.sum(res['track']['kl'] * alpha) + categorical_kl(alpha, pi)

        # 2. compute psi and psi2 (first and second moments of predictions)
        Xb = data['X'] @ (alpha * res['params']['mu'])
        Zd = jnp.sum(alpha * res['params']['delta']) 
        psi = Xb + Zd

        b2 = alpha * (res['params']['mu']**2 + 1/res['params']['tau'])
        Xb2 = data['X']**2 @ b2
        VXb = Xb2 - Xb**2
        psi2 = psi**2 + VXb

        # 3. compute likelihood under null model
        null_model = tilted_null_model(data['y'], re)
        null_likelihood = null_model['track']['elbo']

        # 4. compute lbf per feature and whole model
        lbf_variable = res['track']['elbo'] - null_likelihood
        lbf = jnp.sum(lbf_variable * alpha) - categorical_kl(alpha, pi) 

        # package up parameters and summary of SER
        res['params']['alpha'] = alpha
        res['params']['pi'] = pi

        # add ser elbo to track
        res['track']['elbo_ser'] = elbo
        res['track']['kl_ser'] = kl

        summary = dict(
            lbf_variable = lbf_variable,
            lbf = lbf,
            psi = psi,
            psi2 = psi2
        )
        res['summary'] = summary
        return res
    return fit_ser_jax

fit_uvb_ser = make_ser(fit_univariate_vb)
fit_tilted_ser = make_ser(tilted_univariate_lr)

def initialize_ser(X, y, tau0, o, o2):
    X = np.array(X)
    y = np.array(y)
    data = dict(X=X, y=y)
    n, p = X.shape
    
    
    re = initialize_re(n, o, o2)
    params = initialize_ser_params(n, p, tau0)
    return(dict(data=data, re=re, params=params, control={}))

def fit_uvb_ser2(X, y, tau0, o, o2):
    return map_nested_dicts(
      fit_uvb_ser(**initialize_ser(X, y, tau0, o, o2)),
      np.array
    )

def fit_tilted_ser2(X, y, tau0):
    return map_nested_dicts(
      fit_tilted_ser(**initialize_ser(X, y, tau0)),
      np.array
    )
