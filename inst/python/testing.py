from inst.python.ibss2m import ibss2m_uvb
from inst.python.ser import fit_uvb_ser, initialize_re, initialize_ser_params, initialize_univariate_params
from inst.python.univariate_vb import fit_univariate_vb
import numpy as np

def simulate_data(n, p):
    X = np.random.normal(size=n*p).reshape(n, -1)
    y = np.random.binomial(1, 1/(1 + np.exp(1 - 2* X[:, :1].sum(1))), n)
    return dict(X=X, y=y)

def get_data_x(data, idx):
    return dict(X = data['X'][:, idx], y=data['y'])

n = 1000
p = 500

re = initialize_re(n, mu2=1)
params = initialize_ser_params(n, p, 1.0)
data = simulate_data(n, p)

# fit single
params = initialize_univariate_params(n, 1.0)
fit1 = fit_univariate_vb(get_data_x(data, 0), re, params)
fit1['params'].pop('xi')
fit1['params']

# fit SER
re = initialize_re(n, mu2=1)
params = initialize_ser_params(n, p, 1.0)
fit3 = fit_uvb_ser(data, re, params, {})

# fit IBSS2m
ibss_fit = ibss2m_uvb(data['X'], data['y'], L=3)
ibss_fit2 = ibss2m_uvb(data['X'], data['y'], L=3, keep_2m=False, tol=1e-6)

