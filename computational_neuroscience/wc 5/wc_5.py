import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal

N = 10000
sigma_2 = 1
rho = 0.99
mu = [-1,2]
Q = [[1,rho],[rho,1]]


X,Y = np.mgrid[-3.5:1.5:.01, -0.5:4.5:.01]
pos = np.dstack([X,Y])
normal_distr = multivariate_normal(mu,Q)
plt.contourf(X,Y,normal_distr.pdf(pos))





