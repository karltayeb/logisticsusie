import numpy as np
import jax.numpy as jnp
from jax.scipy.special import logsumexp

def sigmoid(x):
    return 0.5 * (jnp.tanh(x / 2) + 1)

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
