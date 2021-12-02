import numpy as np
import matplotlib.pyplot as plt

def dphase(K,phase,freq,R,Theta):
    d_phase = freq-K*R*np.sin(phase-Theta)
    return d_phase 

def mean_field(freq): #calculates the mean field of a bunch of oscillators
    r_exp_iO = np.sum(np.exp(freq*1j))/N
    R = np.abs(r_exp_iO)
    Theta = np.angle(r_exp_iO)
    return R, Theta

def modified_euler(K,phase,freq,R,Theta,h): #Modified Euler algorithm to get new phase 
    phase_tilde = phase + h*dphase(K,phase,freq,R,Theta)
    phase_new = phase + 0.5*h*(dphase(K,phase,freq,R,Theta)+dphase(K,phase_tilde,freq,R,Theta))
    return phase_new

def evolution(K,phase,freq,h=0.1,T=1000):
    phase_list = []
    R_list = []
    Theta_list = []
    for t in range(int(T/h)):
        print(round(t/int(T/h),2))
        R,Theta = mean_field(phase0)
        phase = modified_euler(K,phase,freq,R,Theta,h)
        phase_list.append(phase)
        R_list.append(R)
        Theta_list.append(Theta)
    return np.array(phase_list).transpose(),np.array(R_list).transpose(),np.array(Theta_list).transpose()

N = 1000
freq = np.random.normal(0,0.1,size=N) #oscillator freq
phase0 = np.random.normal(0,1,size=N) #starting phase of oscillator 
K = 1

phase,R,Theta = evolution(K,phase0,freq)
phase = (phase+np.pi)%(2*np.pi)-np.pi
Theta = (Theta+np.pi)%(2*np.pi)-np.pi



#plotting R,Theta for different values of K
K = [[0,0.05,0.1],[0.15,0.2,0.25,0.3]]
N = 10 #adapt this to get different amount of neurons 
T = 20
freq = np.random.normal(0,0.1,size=N) #oscillator freq
phase0 = np.random.normal(0,1,size=N) #starting phase of oscillator 


#phase plots 
fig1,ax1 = plt.subplots(nrows=2,ncols=3)
fig2,ax2 = plt.subplots(nrows=2,ncols=3)
fig1.tight_layout()
fig2.tight_layout()
for i in range(2):
    for j in range(3):
        t,phase,R,Theta = evolution(K[i][j],phase0,freq,T)
        phase = (phase+np.pi)%(2*np.pi)-np.pi
        Theta = (Theta+np.pi)%(2*np.pi)-np.pi
        for m in range(N):
            ax1[i][j].plot(t,phase[m])
        ax1[i][j].plot(t,Theta,label="Big Theta")
        ax2[i][j].plot(t,R)
        if i == 1:
            ax1[i][j].set_xlabel("time")
            ax2[i][j].set_xlabel("time")
        elif i == 0:   
            ax1[i][j].set_xticks([])
            ax2[i][j].set_xticks([])
        if j == 0:
            ax1[i][j].set_ylabel("phase")
            ax2[i][j].set_ylabel("phase")
        ax1[i][j].legend(loc="upper center")
        ax1[i][j].set_title("K = {}, N={}".format(K[i][j],N))
        ax1[i][j].set_ylim([-np.pi,np.pi+2]) #to show the legend the right way 
        ax2[i][j].set_title("K = {}, N={}".format(K[i][j],N))
        ax2[i][j].set_ylim([0,1.3])
