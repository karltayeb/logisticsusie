from inst.python.ibss2m import add_re, ibss2m_tilted
from inst.python.ser import fit_tilted_ser, initialize_re, initialize_ser_params
import numpy as np

def simulate_data(n, p):
    X = np.random.normal(size=n*p).reshape(n, -1)
    y = np.random.binomial(1, 1/(1 + np.exp(1 - 2* X[:, :1].sum(1))), n)
    return dict(X=X, y=y)

def get_data_x(data, idx):
    return dict(X = data['X'][:, idx], y=data['y'])

# simulate data
n = 1000
p = 500
re = initialize_re(n, mu2=1)
params = initialize_ser_params(n, p, 1.0)
data = simulate_data(n, p)

# fit ser
params = initialize_ser_params(n, p, 10)
re = initialize_re(n, mu2=1)
ser = fit_tilted_ser(data, re, params, {})
ser['summary']['lbf']

# fit ibss
ibss_fit = ibss2m_tilted(data['X'], data['y'], L=3)
ibss_fit = ibss2m_tilted(data['X'], data['y'], L=3, keep_2m=False, tol=1e-6)