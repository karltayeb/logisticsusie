import numpy as np
import jax.numpy as jnp
from jax.scipy.special import logsumexp
import matplotlib.pyplot as plt
import collections

def map_nested_dicts(ob, func):
    """
    Apply func to each element of a nested dict, recursively
    """
    if isinstance(ob, dict):
        return {k: map_nested_dicts(v, func) for k, v in ob.items()}
    else:
        return func(ob)

def sigmoid(x):
    return 0.5 * (jnp.tanh(x / 2) + 1)

def logit(x):
    return 2 * jnp.arctanh(2*x - 1)

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

def plot_contour(f, xrange: tuple = (-10, 10, 50), yrange: tuple = (-10, 10, 50)):
    x = np.linspace(*xrange)
    y = np.linspace(*yrange)

    X, Y = np.meshgrid(x, y)
    Z = f(X, Y)

    plt.contourf(X, Y, Z, 20, cmap='RdGy')
    plt.colorbar()
    plt.show()
    plt.close()

def plot_1d(f, xrange: tuple= (-10, 10, 50)):
    x = np.linspace(*xrange)
    y = f(x)
    plt.plot(x, y)
    plt.show()
    plt.close()
