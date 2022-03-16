import numpy as np
import matplotlib.pyplot as plt

# Bellman equations to calculate all backward values 

def create_states(h):
    N = int((h**4+10*h**3+35*h**2+50*h+24)/24) # total amount of possible states evaluated from 4 double sums 
    all_states = np.zeros([N,5],dtype=int) # empty array of len V_len
    count = 0
    for t in range(h+1): # look at all time points 
        for n1 in range(t+1): # look at all states at a certain time
            n2 = t-n1 # (n1,n2) = (10,0),(9,1),(8,2),...,(9,0),(8,1),(7,2),...,(8,0),(7,1),...,...(0,1),(1,0),(0,0) --> state (0,0) can maybe cause problems       
            for w1 in range(n1+1):
                for w2 in range(n2+1):
                    all_states[count] = t,n1,w1,n2,w2
                    count += 1
    return all_states

def get_index(n1,w1,n2,w2,all_states): # get index that correlates with certain state from the 'all_state' matrix
    if type(n1)==int or type(n1)==float:
        index_1 = np.intersect1d(np.where(all_states[:,1]==n1),np.where(all_states[:,2]==w1))
        index_2 = np.intersect1d(np.where(all_states[:,3]==n2),np.where(all_states[:,4]==w2))
        index = np.intersect1d(index_1,index_2)[0]
    else:
        index = np.zeros(len(n1),dtype=int)
        for i in range(len(n1)):
            index_1 = np.intersect1d(np.where(all_states[:,1]==n1[i]),np.where(all_states[:,2]==w1[i]))
            index_2 = np.intersect1d(np.where(all_states[:,3]==n2[i]),np.where(all_states[:,4]==w2[i]))
            index_i = np.intersect1d(index_1,index_2)[0]
            index[i] = index_i
    return index 

def p(n,w): # calculate expected probability from drawings and payoffs for a machine
    return (w+1)/(n+2)

def get_V(all_states):
    N = all_states.shape[0]
    V = np.zeros(N)

    index = np.where(all_states[:,0]==h-1)
    t_states = all_states[index]
    n1,w1,n2,w2 = t_states[:,1],t_states[:,2],t_states[:,3],t_states[:,4]
    p1,p2 = p(n1,w1),p(n2,w2)
    V[index] = np.max([p1,p2],axis=0)

    for t in range(h-1)[::-1]: 
        index = np.where(all_states[:,0]==t)
        t_states = all_states[index]
        n1,w1,n2,w2 = t_states[:,1],t_states[:,2],t_states[:,3],t_states[:,4]
        p1,p2 = p(n1,w1),p(n2,w2)
        value1 = p1*(1+V[get_index(n1+1,w1+1,n2,w2,all_states)])+(1-p1)*V[get_index(n1+1,w1,n2,w2,all_states)]
        value2 = p2*(1+V[get_index(n1,w1,n2+1,w2+1,all_states)])+(1-p2)*V[get_index(n1,w1,n2+1,w2,all_states)]
        V[index] = np.max([value1,value2],axis=0)
    return V

h = 4
all_states = create_states(h)
V = get_V(all_states)
V[get_index(0,0,1,0,all_states)] # This should be 1.53 as stated in the slides for the state (0010)



#####################################

# initialize values before drawing
p1,p2 = 0.3,0.5 # probability of machine 1,2 for payoff 1
n1,w1,n2,w2 = 0,0,0,0 # amount of draws and payoffs for machine 1,2 in the beginning
h = 10 # total amount of draws
state = np.zeros([h,4]) # state matrix
n1 += 1 # draw from machine 1 in first step
w1 = 1 if np.random.rand() <= p1 else 0 # payoff with probability 1
state[1] = [n1,w1,n2,w2] # update first step in the state
# Start drawing other samples using all possible states V*
# 1) Calculate expected outcome after the first step via argmax and pick the highest probability
# 2) Calculate expected outcome after the second step via argmax and pick the highest probability
# etc...
# In the end you have a matrix states with in the end a state: states[h] that gives the final outcome. This immediately gives the final payoff w=w1+w2