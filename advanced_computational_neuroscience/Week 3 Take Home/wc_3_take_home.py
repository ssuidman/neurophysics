import numpy as np
import matplotlib.pyplot as plt

# data = np.loadtxt("/Users/samsuidman/Desktop/neurophysics/advanced_computational_neuroscience/Week 3 Take Home/salamander_data.txt")

s_all = data.reshape(160,297,953) # s_all[neurons,experiments,time] gives for all neurons all experiments with all timestamps
random = False

if random:    
    neurons = np.sort(np.random.randint(0,160,size=10)) # 10 random neurons sorted on index
    while(len(np.unique(neurons))!=10): # check if there are not double indices
        neurons = np.sort(np.random.randint(0,160,size=10)) # indices of 10 random neurons out of the total
if not random:
    neurons = np.arange(10) # indices of the first 10 of the total amount neurons

s = s_all[neurons,0,:] # s[neurons,patterns] choose 10 neurons and the zero'th experiment 
N = s.shape[0] # amount of neurons
P = s.shape[1] # amount of patterns (=timestamps) 
s_c = 1/P * np.sum(s,axis=1) # s_c[neurons]     = sum( s[neurons,patterns] ) clamped statistics
ss_c = 1/P * np.matmul(s,np.transpose(s)) # ss_c[neurons,neurons]       = sum( s[neurons,patterns]*s[neurons,patterns]^T ) clamped statistics
w_random = np.random.uniform(-1,1,size=[N,N]) # w_random[neurons,neurons] is random connection matrix  
w_symmetric = w_random + np.transpose(w_random) # w_symmetric[neurons,neurons] is random symmetric connection matrix
w = w_symmetric - np.diag(np.diag(w_symmetric)) # w[neurons,neurons] is random symmetric matrix with zeros on diagonal
theta = np.random.uniform(-1,1,size=N) # theta[neurons]

s_every_pattern = np.zeros([N,2**N]) # s_every_pattern[neurons,patterns]
for i in range(2**N):
    s_every_pattern[:,i] = np.array(list(np.binary_repr(i,width=N))).astype(int) # s_every_pattern[patterns,neurons] stores every possible pattern of N neurons by representing numbers between 0 and 2^N-1 in binary

for i in range(10000):
    s_w_s = np.tensordot(np.tensordot(s_every_pattern,w,axes=(0,0)),s_every_pattern,axes=(1,0)) #s_w_s[patterns,patterns] is a sum over s,w,s for all patterns, you want only the diagonal of this because then you sum over the same pattern.
    theta_s = np.tensordot(theta,s_every_pattern,axes=(0,0)) # theta_s[patterns] sums over all neurons for each pattern
    E = - 0.5*np.diag(s_w_s) - theta_s # calculate E[patterns] for each pattern
    Z = np.sum(np.exp(-E)) # calculate Z by summing over the exponential of E(s) for all patterns
    p = 1/Z * np.exp(-E) # p[patterns] gives the probability of each of the possible patterns of s_every_pattern. Therefore sum(p)=1
    s_i = np.matmul(s_every_pattern,p) # s_i[neurons] free statistic variable
    ss_ij = np.matmul(np.multiply(s_every_pattern,p),np.transpose(s_every_pattern)) # ss_ij[neurons,neurons] free statistic variable

    w += learning_rate * (ss_c-ss_ij)
    theta += learning_rate * (s_c-s_i)
    if i in 1000*np.arange(10):
        print(np.sum(np.abs(s_c-s_i)))
        print(np.sum(np.abs(ss_c-ss_ij)))
        print('\n')