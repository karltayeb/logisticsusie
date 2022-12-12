import jax.numpy as jnp
import jax
from jax import vmap, jit
from jax.lax import while_loop
import numpy as np
from jax.scipy.special import logsumexp

def polya_gamma_mean(b, c):
    """
    Mean of PG(b, c)
    """
    # TODO: PG(b, 0) = b/4, deal with divide by 0
    # return b/c * (sigmoid(c) - 0.5)  # equivalent 
    return 0.5 * b/c * jnp.tanh(c/2)

def normal_kl(mu, var, mu0=0, var0=1):
    return 0.5 * (jnp.log(var0) - jnp.log(var) + var/var0 + (mu - mu0)**2/var0 - 1)

def categorical_kl(alpha, pi):
    kl = jnp.nansum(alpha * (jnp.log(alpha) - jnp.log(pi)))
    return kl

def update_intercept(x, y, mu, tau, xi, delta, tau0, offset, offset2):
    kappa = y - 0.5
    xb = x * mu + offset
    omega = polya_gamma_mean(1, xi)
    delta = jnp.sum(kappa - xb * omega) / jnp.sum(omega)
    return delta

def update_b(x, y, mu, tau, xi, delta, tau0, offset, offset2):
    omega = polya_gamma_mean(1, xi)
    kappa = y - 0.5
    tau = jnp.sum(omega * x**2) + tau0
    nu = jnp.sum((kappa - omega * (delta + offset)) * x)
    return dict(mu = nu/tau, tau=tau)

def compute_psi(x, y, mu, tau, xi, delta, tau0, offset, offset2):
    xb = x * mu
    psi = xb + delta + offset
    return psi
  
def compute_psi2(x, y, mu, tau, xi, delta, tau0, offset, offset2):
    Vo = offset2 - offset**2
    Vxb = x**2 / tau
    xb = x * mu
    psi2 = (xb + delta + offset)**2 + Vo + Vxb
    return psi2

def update_xi(x, y, mu, tau, xi, delta, tau0, offset, offset2):
    psi2 = compute_psi2(x, y, mu, tau, xi, delta, tau0, offset, offset2)
    xi = jnp.sqrt(psi2)
    return xi

def compute_elbo(x, y, mu, tau, xi, delta, tau0, offset, offset2):
    kappa = y - 0.5
    psi = compute_psi(x, y, mu, tau, xi, delta, tau0, offset, offset2)
    bound = jnp.log(jax.nn.sigmoid(xi)) + (kappa * psi) - 0.5 * xi
    kl = normal_kl(mu, 1/tau, 0, 1/tau0)
    
    elbo = jnp.sum(bound) - kl
    return elbo

def vb_iter(val):
    # possible improved organization
    # data = x, y, offset, offset2
    # variational_params = mu, tau, xi
    # other = delta, tau0, elbo, diff, it
    
    # unpack
    data = val['data']
    params = val['params']
    re = val['re']
    track = val['track']

    x = data['X']
    y = data['y']

    mu = params['mu']
    tau = params['tau']
    xi = params['xi']
    delta = params['delta']
    tau0 = params['tau0']

    offset = re['mu']
    offset2 = re['mu2']

    # do updates
    delta = update_intercept(x, y, mu, tau, xi, delta, tau0, offset, offset2)
    
    b = update_b(x, y, mu, tau, xi, delta, tau0, offset, offset2)
    mu = b['mu']
    tau = b['tau']
    
    xi = update_xi(x, y, mu, tau, xi, delta, tau0, offset, offset2)
    elbo_new = compute_elbo(x, y, mu, tau, xi, delta, tau0, offset, offset2)
    
    diff = elbo_new - track['elbo']
    it = track['it'] + 1
    # val = (x, y, mu, tau, xi, delta, tau0, offset, offset2, elbo_new, diff, it)

    new_params = dict(mu=mu, tau=tau, xi=xi, delta=delta, tau0=tau0)
    new_track = dict(diff = diff, elbo = elbo_new, it = it)
    new_val = dict(data = data, re=re, params = new_params, track=new_track)
    return new_val

@jit
def fit_univariate_vb(data: dict, re: dict, params: dict):
    def vb_cond(val):
        return val['track']['diff'] > 1e-6
    # initialize
    # data = dict(x=x, y=y)
    # re = dict(mu=offset, mu2=offset2)
    # params = dict(mu=mu, tau=tau, xi=xi, delta=delta, tau0=tau0)
    track = dict(elbo = -jnp.inf, diff=jnp.inf, it=0)
    val_init = dict(data=data, re=re, params=params, track=track)
    
    # iterate updates until convergence
    res = while_loop(vb_cond, vb_iter, val_init)
    # return results

    out = dict(
        params = res['params'],
        track = res['track']
    )
    return out 

# vectorize fit_univariatet VB over X
# Note that othere data (y, offset, and offset2) are not mapped over
univariate_vb_vec_jax = jit(vmap( 
    fit_univariate_vb, 
    in_axes=({'X': 1, 'y': None}, {'mu': None, 'mu2':None}, {'mu': 0, 'tau': 0, 'xi': 1, 'delta': 0, 'tau0': 0}),
    out_axes={'params': {'mu': 0, 'tau': 0, 'xi': 1, 'delta': 0, 'tau0': 0}, 'track': {'elbo': 0, 'diff': 0, 'it': 0}}
))

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


def fit_null_model(y: jnp.array, re: dict):
    """Fit null model where the effect is exactly 0

    Args:
        y (jnp.array): n vector of binary observations
        re (dict): random effects term a 'mu' and 'mu2'
    """
    n = y.size
    X = y[::-1]
    data = dict(X=X, y=y)
    params = initialize_univariate_params(n, 1e10)
    null_model = fit_univariate_vb(data, re, params)
    return null_model

@jit
def fit_uvb_ser_jax(data: dict, re: dict, params: dict, control : dict):
    """
    data: dict with X and y
    re: dict with mu and mu2, first and second moments of random effect per observation
    params: dict with fitting parameters
    control: dict estimate_prior_variance
    """
    # 0: compute UVB for each observation
    pi = params.pop('pi')  # the univariate VBs don't know about pi
    alpha = params.pop('alpha')
    res = univariate_vb_vec_jax(data, re, params)

    # 1. compute pips (alpha) by normalizing elbos
    logits = res['track']['elbo'] + jnp.log(pi)  # unnormalized
    alpha = jnp.exp(logits - logsumexp(logits))
    elbo = np.sum(res['track']['elbo'] * alpha) - categorical_kl(alpha, pi)

    # 2. compute psi and psi2 (first and second moments of predictions)
    Xb = data['X'] @ (alpha * res['params']['mu'])
    Zd = jnp.sum(alpha * res['params']['delta']) 
    psi = Xb + Zd

    b2 = alpha * (res['params']['mu']**2 + 1/res['params']['tau'])
    Xb2 = data['X']**2 @ b2
    VXb = Xb2 - Xb**2
    psi2 = psi**2 + VXb

    # 3. compute likelihood under null model
    null_model = fit_null_model(data['y'], re)
    null_likelihood = null_model['track']['elbo']

    # 4. compute lbf per feature and whole model
    lbf_variable = res['track']['elbo'] - null_likelihood
    lbf = jnp.sum(lbf_variable * alpha) - categorical_kl(alpha, pi) 

    # package up parameters and summary of SER
    res['params']['alpha'] = alpha
    res['params']['pi'] = pi

    summary = dict(
        lbf_variable = lbf_variable,
        lbf = lbf,
        psi = psi,
        psi2 = psi2
    )
    res['summary'] = summary
    return res


def fit_uvb_ser_jax2(X, y, o = None, prior_variance=1.0):
    n, p = X.shape
    if o is None:
        re = initialize_re(n)
    else:
        re = o

    params = initialize_ser_params(n, p, 1/prior_variance)
    data = dict(X = np.array(X), y=np.array(y))
    fit = fit_uvb_ser_jax(data, re, params, {})

    res = {k: np.array(v) for k, v in fit['params'].items()}
    res.pop('xi')
    res.update({k: np.array(v) for k, v in fit['summary'].items()})
    res.update({k: np.array(v) for k, v in fit['track'].items()})
    return res