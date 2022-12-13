from univariate_vb import initialize_re, initialize_ser_params, fit_univariate_vb, univariate_vb_vec_jax, fit_uvb_ser_jax, initialize_univariate_params
from ibss2m import ibss2m_jax
import numpy as np

def simulate_data(n, p):
    X = np.random.normal(size=n*p).reshape(n, -1)
    y = np.random.binomial(1, 1/(1 + np.exp(-X[:, :2].sum(1))), n)
    return dict(X=X, y=y)

def get_data_x(data, idx):
    return dict(X = data['X'][:, 0], y=data['y'])

n = 1000
p = 500

re = initialize_re(n, mu2=1)
params = initialize_ser_params(n, p, 1.0)
data = simulate_data(n, p)

# fit single
params = initialize_univariate_params(n, 1.0)
fit1 = fit_univariate_vb(get_data_x(data, 1), re, params)

# fit vectorized
re = initialize_re(n, mu2=1)
params = initialize_ser_params(n, p, 1.0)
params.pop('pi')
params.pop('alpha')
fit2 = univariate_vb_vec_jax(data, re, params)

# fit SER
re = initialize_re(n, mu2=1)
params = initialize_ser_params(n, p, 1.0)
fit3 = fit_uvb_ser_jax(data, re, params, {})

# fit IBSS2m
ibss_fit = ibss2m_jax(data['X'], data['y'], L=3)
