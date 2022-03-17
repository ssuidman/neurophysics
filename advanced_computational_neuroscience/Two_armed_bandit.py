import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Bellman equations to calculate all backward values 

def create_states(h):
    N = int((h**4+10*h**3+35*h**2+50*h+24)/24) # total amount of possible states evaluated from 4 double sums 
    all_states = np.zeros([N,5],dtype=int) # empty array of len V_len
    index_df = {}
    count = 0
    for t in range(h+1): # look at all time points 
        for n1 in range(t+1): # look at all states at a certain time
            n2 = t-n1 # (n1,n2) = (10,0),(9,1),(8,2),...,(9,0),(8,1),(7,2),...,(8,0),(7,1),...,...(0,1),(1,0),(0,0) --> state (0,0) can maybe cause problems       
            for w1 in range(n1+1):
                for w2 in range(n2+1):
                    all_states[count] = t,n1,w1,n2,w2
                    index_df[(n1,w1,n2,w2)] = count
                    count += 1
    index_df = pd.DataFrame([index_df])
    return all_states,index_df

def get_index(index_array,index_df):
    index_array = np.array(index_array)
    if index_array.ndim==1:
        i = index_df[tuple(index_array)].values[0]
    else:
        i = index_df[pd.Index(tuple(map(tuple,index_array)))].values[0]
    return i

def p(n,w): # calculate expected probability from drawings and payoffs for a machine
    return (w+1)/(n+2)

def get_V(all_states,index_df):
    N = all_states.shape[0]
    V = np.zeros(N)

    i = np.where(all_states[:,0]==h-1) # do "last" step before the rest, because V is calculated from p1,p2 only
    t_states = all_states[i][:,1:]
    n1,w1,n2,w2 = t_states[:,0],t_states[:,1],t_states[:,2],t_states[:,3]
    p1,p2 = p(n1,w1),p(n2,w2)
    V[i] = np.max([p1,p2],axis=0)

    for t in range(h-1)[::-1]: 
        i = np.where(all_states[:,0]==t)
        t_states = all_states[i][:,1:]
        n1,w1,n2,w2 = t_states[:,0],t_states[:,1],t_states[:,2],t_states[:,3]
        V1 = get_index(t_states+np.array([1,1,0,0]),index_df)
        V2 = get_index(t_states+np.array([1,0,0,0]),index_df)
        V3 = get_index(t_states+np.array([0,0,1,1]),index_df)
        V4 = get_index(t_states+np.array([0,0,1,0]),index_df)
        value1 = p(n1,w1)*(1+V[V1])+(1-p(n1,w1))*V[V2]
        value2 = p(n2,w2)*(1+V[V3])+(1-p(n2,w2))*V[V4]
        V[i] = np.max([value1,value2],axis=0)
    return V

h = 4
all_states,index_df = create_states(h)
V = get_V(all_states,index_df)
V[index_df[(1,0,1,0)].values[0]]


#################################################################################################
##################################### Start the random walk #####################################
#################################################################################################


h = 10
all_states,index_df = create_states(h)
V = get_V(all_states,index_df)
# initialize values before drawing
p1,p2 = 0.5,0.5 # probability of machine 1,2 for payoff 1
n1,w1,n2,w2 = 0,0,0,0 # initialize states
states = np.zeros([h+1,4]) # state matrix

for i in range(1,h+1): # first state is all zeros, from there fill in the matrix
    V1 = get_index([n1+1,w1+1,n2,w2],index_df)
    V2 = get_index([n1+1,w1,n2,w2],index_df)
    V3 = get_index([n1,w1,n2+1,w2+1],index_df)
    V4 = get_index([n1,w1,n2+1,w2],index_df)
    value1 = p(n1,w1)*(1+V[V1])+(1-p(n1,w1))*V[V2]
    value2 = p(n2,w2)*(1+V[V3])+(1-p(n2,w2))*V[V4]
    machine =  np.argmax([value1,value2]) + 1
    if machine==1:
        n1 += 1
        if np.random.rand()<p1:
            w1 += 1
    else:
        n2 += 1
        if np.random.rand()<p2:
            w2 += 1
    states[i] = n1,w1,n2,w2
states
