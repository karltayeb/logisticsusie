import jax.numpy as jnp
import numpy as np
from jax import grad, hessian
from jax.lax import while_loop
from jax import jit, vmap
from jax.scipy.special import logsumexp
from inst.python.utils import sigmoid, normal_kl

def tilted_psi(data: dict, re: dict, params: dict):
    psi0 = re['mu']
    Vpsi0 = re['mu2'] - re['mu']**2

    psi = params['delta'] + data['X'] * params['mu']
    Vpsi = data['X']**2 * params['var']

    return psi + psi0, Vpsi + Vpsi0

def tilted_bound(data: dict, re: dict, params: dict, params2: dict):
    psi, Vpsi = tilted_psi(data, re, params)
    tmp = -(psi + 0.5 * (1 - 2 * params2['xi']) * Vpsi)
    bound = data['y'] * psi \
        - 0.5 * params2['xi']**2 * Vpsi \
        + jnp.log(sigmoid(tmp))
    bound = jnp.sum(bound)
    return bound

def xi_fixed_point(data: dict, re: dict, params: dict, params2: dict):
    psi, Vpsi = tilted_psi(data, re, params)
    xi = params2['xi']
    for _ in range(10):
        tmp = psi + 0.5 * (1 - 2 * xi) * Vpsi
        xi = sigmoid(tmp)
    return xi

def tilted_elbo(data: dict, re: dict, params: dict, params2: dict):
    elbo = tilted_bound(data, re, params, params2) \
        - normal_kl(params['mu'], params['var'], 0, 1/params2['tau0'])
    return elbo

# convert parameter dictionary to unconstrained parameter vector
param2dict = lambda x: dict(mu = x[0], var=jnp.exp(x[1]), delta=x[2])
dict2param = lambda d: jnp.array([d['mu'], jnp.log(d['var']), d['delta']])

# objective function with mu, var, and delta in vector format
def tilted_obj(x, data, re, params2):
    params = param2dict(x)
    return tilted_elbo(data, re, params, params2)

# gradient and hessian for newton method
tilted_grad = grad(tilted_obj)
tilted_H = hessian(tilted_obj)

def cond_fun(val):
    # not_converged = jnp.logical_and(
    #     jnp.abs(val['diff']) > val['tol'],
    #     val['it'] < val['maxit']
    # )
    not_converged = jnp.abs(val['diff']) > val['tol']
    return not_converged

def body_fun(val):
    params = val['params']
    params2 = val['params2']
    data = val['data']
    re = val['re']
    old_elbo = val['elbo']
    tol = val['tol']
    it = val['it']
    maxit = val['maxit']

    gradStep = lambda x: x - jnp.linalg.solve(tilted_H(
        x, data, re, params2), tilted_grad(x, data, re, params2))
    params = gradStep(params)
    params2['xi'] = xi_fixed_point(data, re, param2dict(params), params2)

    elbo = tilted_obj(params, data, re, params2)
    diff = elbo - old_elbo

    val = dict(
        params = params,
        params2 = params2,
        data = data,
        re = re,
        diff = diff,
        elbo = elbo,
        tol = tol,
        it = it + 1,
        maxit = maxit
    )
    return val

def tilted_univariate_lr(data, re, params):
    # break up params into optimization groups
    params2 = params
    params = dict(
        mu = params2.pop('mu'),
        var = 1/params2.pop('tau'),
        delta = params2.pop('delta'),
    )
    params = dict2param(params)

    # initialize
    val_init = dict(
        params = params,
        params2 = params2,
        data = data,
        re = re,
        diff = 1000,
        elbo = tilted_obj(params, data, re, params2),
        tol = 1e-4,
        it = 0,
        maxit = 50
    )

    # fit
    val = while_loop(cond_fun, body_fun, val_init)

    # wrapup-- outuput consistent with univariate_vb
    # {
    # 'params': {'mu', 'tau', 'xi', 'delta', 'tau0'},
    # 'track': {'elbo', 'kl', 'diff', 'it'}
    # }

    val['params'] = param2dict(val['params'])
    params = dict(
        mu = val['params']['mu'],
        tau = 1/val['params']['var'],
        xi = val['params2']['xi'],
        delta = val['params']['delta'],
        tau0 = val['params2']['tau0']
    )
    track = dict(
        elbo = val['elbo'],
        kl = normal_kl(params['mu'], 1/params['tau'], 0, 1/params2['tau0']),
        diff = val['diff'],
        it = val['it'] 
    )
    res = dict(params=params, track=track)
    return res
