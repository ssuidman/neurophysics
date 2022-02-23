import numpy as np
import matplotlib.pyplot as plt

# data = np.loadtxt("/Users/samsuidman/Desktop/neurophysics/advanced_computational_neuroscience/Week 3 Take Home/salamander_data.txt")

s_all = data.reshape(160,297,953) # s_all[neurons,experiments,time]
random = False

if random:    
    neurons = np.sort(np.random.randint(0,160,size=10))
    while(len(np.unique(neurons))!=10):
        neurons = np.sort(np.random.randint(0,160,size=10))
if not random:
    neurons = np.arange(10) 

s = s_all[neurons,0,:] # s[neurons,time] choose 10 neurons and the zero'th experiment 
P = s.shape[0]
T = s.shape[1]
s_c = 1/P * np.sum(s,axis=0) # s_c[time]     = sum( s[neurons,time] ) clamped statistics
ss_c = 1/P * np.matmul(np.transpose(s),s) # ss_c[time,time]       = sum( s[neurons,time]^T*s[neurons,time] ) clamped statistics
w = np.random.uniform(-1,1,size=[T,T]) # w[time,time]
theta = np.random.uniform(-1,1,size=T) # theta[time]


E = -0.5*np.diag(np.matmul(np.matmul(s,w),np.transpose(s))) - np.matmul(s,theta) # here things start to go wrong 
p = np.exp(-E) # p[neurons] 
# Z = something to normalize p
# s_i =  # s_i[time]  free statistics
# ss_ij # ss_ij[time,time]  free statistics
