import numpy as np
import matplotlib.pyplot as plt
from time import time

def p(n,w): # calculate expected probability from drawings and payoffs for a machine
    return (w+1)/(n+2)

def get_V(h,p1,p2):
    t_array = np.zeros(h)
    V = np.zeros([h+1,h+1,h+1,h+1])
    print('h =',h,', p1 =',p1,', p2 =',p2)
    print('create V with N =',int((h**4+10*h**3+35*h**2+50*h+24)/24),'states')
    t1 = time()
    count = 0
    for n1 in range(h): # look at all states at a certain time
        n2 = h-1-n1 # (n1,n2) = (10,0),(9,1),(8,2),...,(9,0),(8,1),(7,2),...,(8,0),(7,1),...,...(0,1),(1,0),(0,0) --> state (0,0) can maybe cause problems       
        for w1 in range(n1+1):
            for w2 in range(n2+1):
                count += 1
                V[n1,w1,n2,w2] = np.max([p(n1,w1),p(n2,w2)])
    t2 = time()
    print('t =',h,round(t2-t1,4),'sec for',count,'states')
    t_array[h-1] = t2-t1

    for t in range(h)[::-1]: # look at all time points 
        count = 0
        t1 = time()
        for n1 in range(t+1): # look at all states at a certain time
            n2 = t-1-n1 # (n1,n2) = (10,0),(9,1),(8,2),...,(9,0),(8,1),(7,2),...,(8,0),(7,1),...,...(0,1),(1,0),(0,0) --> state (0,0) can maybe cause problems       
            for w1 in range(n1+1):
                for w2 in range(n2+1):
                    count += 1
                    value1 = p(n1,w1)*(1+V[n1+1,w1+1,n2,w2])+(1-p(n1,w1))*V[n1+1,w1,n2,w2]
                    value2 = p(n2,w2)*(1+V[n1,w1,n2+1,w2+1])+(1-p(n2,w2))*V[n1,w1,n2+1,w2]
                    V[n1,w1,n2,w2] = np.max([value1,value2],axis=0)
        t2 = time()
        print('t =',t,round(t2-t1,4),'sec for',count,'states')
        t_array[h-1-t] = t2-t1
    return V,t_array

def two_armed_bandit(h,p1,p2):
    V,t = get_V(h,p1,p2)

    # initialize values before drawing
    n1,w1,n2,w2 = 0,0,0,0 # initialize states
    states = np.zeros([h+1,4]) # state matrix

    for i in range(1,h+1): # first state is all zeros, from there fill in the matrix
        value1 = p(n1,w1)*(1+V[n1+1,w1+1,n2,w2])+(1-p(n1,w1))*V[n1+1,w1,n2,w2]
        value2 = p(n2,w2)*(1+V[n1,w1,n2+1,w2+1])+(1-p(n2,w2))*V[n1,w1,n2+1,w2]
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
        w = w1+w2
    return states,w,V,t

h = 4
p1,p2 = 0.2,0.3 # probability that machine pays off 1 

states,w,V,t = two_armed_bandit(h,p1,p2)
print('Machines chosen:')
print(np.argmax(np.diff(states[:,[0,2]],axis=0),axis=1)+1)
print('Total profit:',w)

if h>20:
    fig,ax = plt.subplots()
    ax.plot(np.arange(h)[::-1],t)
    ax.set_xlabel('t')
    ax.set_ylabel('time per iteration')
    ax.set_xticks(np.arange(0,h,int(h/10)))

print('\n')
print('Example: V =',round(V[1,0,1,0],2),'for h=4')
