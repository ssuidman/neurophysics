import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from scipy.sparse.csgraph import dijkstra
import itertools

# problem 2a
def problem_2a(n,m,printen,plotten):
    A = np.zeros(n**2)
    A[:m] = 1
    np.random.shuffle(A)
    A = np.reshape(A,[n,n])
    A_sum = np.sum(A,axis=0)
    C = A_sum/n
    if printen:
        print(A)
    if plotten:
        plt.figure()
        plt.hist(C,bins=10)
        plt.show()
    return A

def problem_2b(n,printen,plotten):
    A = np.random.choice([0,1],size=[n,n],p=[770/870,100/870])
    A_sum = np.sum(A,axis=0)
    C = A_sum/n
    if printen:
        print(A)
    if plotten:
        plt.figure()
        plt.hist(C,bins=10)
        plt.show()
        return A

def problem_2c(n,k,printen):
    A = np.array([[1 if i!=j and (np.abs(i-j)<=int(k/2) or np.abs(i-j)>=np.abs(n-int(k/2))) else 0 for i in range(n)] for j in trange(n)])
    if printen:
        print(A)
    return A

def problem_2d(A,p,printen):    
    indices_1 = np.array(np.where(A==1)) # indices where the matrix is 1
    n_shuffle = int(indices_1.shape[1]*p/2) # rewire a fraction of these indices, because the matrix is symmetrical there are need to be shuffled p/2 values for the upple triangle
    upper_indices_1 = np.array([indices_1[0,np.where(indices_1[1]>indices_1[0])][0],indices_1[1,np.where(indices_1[1]>indices_1[0])][0]]) # look only at upper triangle indices
    n_1 = upper_indices_1.shape[1] #look at amount of upper triangle indices
    i_1 = np.random.choice(np.arange(n_1),size=n_shuffle,replace=False) # pick n_shuffle random indices where the matrix is 1
    x_1,y_1 = upper_indices_1[0,i_1],upper_indices_1[1,i_1] # get the x,y-index of the elements to be shuffled

    indices_0 = np.array(np.where(A==0)) # indices where the matrix is 1
    upper_indices_0 = np.array([indices_0[0,np.where(indices_0[1]>indices_0[0])][0],indices_0[1,np.where(indices_0[1]>indices_0[0])][0]]) # look only at upper triangle indices
    n_0 = upper_indices_0.shape[1] #look at amount of upper triangle indices
    i_0 = np.random.choice(np.arange(n_0),size=n_shuffle,replace=False) # pick n_shuffle random indices where the matrix is 1
    x_0,y_0 = upper_indices_0[0,i_0],upper_indices_0[1,i_0] # get the x,y-index of the elements to be shuffled

    A[x_1,y_1] = 0 # replace ones by zeros
    A[y_1,x_1] = 0 # make the matrix symmetrical
    A[x_0,y_0] = 1 # replace zeros by ones
    A[y_0,x_0] = 1 # make the matrix symmetrical
    if printen:
        print(A)
    return A

def problem_2e(A): 
    for i in trange(n): # look at each node
        neighbours = np.where(A[i]==1)[0] # get all neighbours
        neighbours_combinations = np.array(list(itertools.combinations(neighbours,2))) # get all combinations of neighbours
        if len(neighbours_combinations)!=0: # look if there are more than 2 neighbours
            triangles = A[neighbours_combinations[:,0],neighbours_combinations[:,1]] # get the values of the combinations of neighbours where a 1 means a connection (=triangle) and a 0 not
            C_node = np.sum(triangles)/neighbours_combinations.shape[0] # get the local cluster index
        else: # this is the case if there are less than 2 neighbours
            C_node = 0 # set C_node to 0 if it has only 1 neighbours
    C = np.mean(C_node) # get the average local cluster index
    d_matrix = dijkstra(A) # make shortest path matrix via Dijkstra algorithm
    d = np.mean(d_matrix) # get the mean path length
    return C,d

A1 = problem_2a(n=30,m=100,printen=False,plotten=False)
A2 = problem_2b(n=30,printen=False,plotten=False)
A3 = problem_2c(n=14,k=4,printen=False)
A4 = problem_2d(A=A3.copy(),p=0.2,printen=False)
C,d = problem_2e(A=A4.copy())


############## OPDRACHT 3 ###################
# Define "real" adjacency matrix that we're trying to estimate
n = 100
mu = np.zeros(n)
a = np.random.rand(n) # arbitrary example
sigma = np.outer(a,a) # "real" covariance matrix
h = 0.2 # threshold
A = np.where(sigma>h,1,0)

# Generate samples from multivariate normal distribution
n_samples = 10 
n_trials = 100 
err = 0 
for trial in range(n_trials):
    x = np.random.multivariate_normal(mu,sigma,n_samples) # Sample from distribution
    sigma_est = np.cov(x,rowvar=False) # Estimate the covariance matrix
    A_est = np.where(sigma_est>h,1,0) # Threshold covariance matrix to estimate adjacency matrix
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

print('Probability of correctly estimating a particular entry in the adjacency matrix: {}'.format(1 - err/n_trials/n**2))