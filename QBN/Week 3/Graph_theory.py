import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange

def problem_2c(n,k,printen):
    A = np.array([[1 if i!=j and (np.abs(i-j)<=k or np.abs(i-j)>=np.abs(n-k)) else 0 for i in range(n)] for j in trange(n)])
    if printen:
        print(A)
    return A

def problem_2d(A,p,printen):    
    n = len(A)
    n_shuffle = int(p*n**2)
    n_1 = np.where(A==1)[0].size
    i_1 = np.random.choice(np.arange(n_1),size=n_shuffle,replace=False)
    x_1,y_1 = np.where(A==1)[0][i_1],np.where(A==1)[1][i_1] # pick random indices to replace ones by zeros
    n_0 = np.where(A==0)[0].size
    i_0 = np.random.choice(np.arange(n_0),size=n_shuffle,replace=False)
    x_0,y_0 = np.where(A==0)[0][i_0],np.where(A==0)[1][i_0] # pick random indices to replace zeros by ones
    A[x_1,y_1] = 0 # replace ones by zeros
    A[x_0,y_0] = 1 # replace zeros by ones
    if printen:
        print(A)
    return A

A3 = problem_2c(n=12,k=3,printen=False)
A4 = problem_2d(A=A3.copy(),p=0.05,printen=False)
print(A3-A4)





################## MATLAB CODE VERTAALD ######################

# QBN graph theory: estimate adjacency matrix based on sample covariance matrix

# Define "real" adjacency matrix that we're trying to estimate
n_nodes = 3
mu = np.zeros(n_nodes)
a = np.array([0.5,1.0,0.3]) # arbitrary example
sigma = np.outer(a,a) # "real" covariance matrix
h = 0.2 # threshold
A = np.where(sigma>h,1,0)

# Generate samples from multivariate normal distribution
n_samples = 10
n_trials = 100
err = 0
for trial in range(n_trials):
    # Sample from distribution
    x = np.random.multivariate_normal(mu,sigma,n_samples)
    
    # Estimate the covariance matrix
    sigma_est = np.cov(x,rowvar=False)

    # Threshold covariance matrix to estimate adjacency matrix
    A_est = np.where(sigma_est>h,1,0)
    err += np.sum(np.where(A==A_est,0,1))

# Plot results
plt.figure(1)
plt.imshow(sigma)
plt.colorbar()
plt.title('Covariance matrix')

plt.figure(2)
plt.imshow(A)
plt.colorbar()
plt.title('Adjacency matrix')

plt.figure(3)
plt.imshow(sigma_est)
plt.colorbar()
plt.title('covariance matrix (estimated)')

plt.figure(4)
plt.imshow(A_est)
plt.colorbar()
plt.title('Adjacency matrix (estimated)')
plt.show()

print('Probability of correctly estimating a particular entry in the adjacency matrix: {}\n'.format(1 - err/n_trials/n_nodes**2))
