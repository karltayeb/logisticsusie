import numpy as np
import jax.numpy as jnp
import jax
from jax import vmap, jit
from jax.lax import while_loop
from jax.scipy.special import logsumexp
from utils import polya_gamma_mean, normal_kl, categorical_kl

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
    kl = normal_kl(mu, 1/tau, 0, 1/tau0)
 
    diff = elbo_new - track['elbo']
    it = track['it'] + 1

    new_params = dict(mu=mu, tau=tau, xi=xi, delta=delta, tau0=tau0)
    new_track = dict(diff = diff, elbo = elbo_new, kl=kl, it = it)
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
    track = dict(elbo = -jnp.inf, kl=0, diff=jnp.inf, it=0)
    val_init = dict(data=data, re=re, params=params, track=track)
    
    # iterate updates until convergence
    res = while_loop(vb_cond, vb_iter, val_init)
    # return results

    out = dict(
        params = res['params'],
        track = res['track']
    )
    return out 
