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
        R,Theta = mean_field(phase)
        phase = modified_euler(K,phase,freq,R,Theta,h)
        phase_list.append(phase)
        R_list.append(R)
        Theta_list.append(Theta)
    return np.linspace(0,T,int(T/h)),np.array(phase_list).transpose(),np.array(R_list).transpose(),np.array(Theta_list).transpose()

#N = 1000
#freq = np.random.normal(0,0.1,size=N) #oscillator freq
#phase0 = np.random.normal(0,1,size=N) #starting phase of oscillator 
#K = 1

#t,phase,R,Theta = evolution(K,phase0,freq)
#phase = (phase+np.pi)%(2*np.pi)-np.pi
#Theta = (Theta+np.pi)%(2*np.pi)-np.pi



# Step 1: 
# 
# plotting R,Theta for different values of K
K = [[0,0.05,0.1],[0.15,0.2,0.25]]
N = 1000 #adapt this to get different amount of neurons 
T = 20 #Take T=100 for Step 2 to see that T=50 most of the time r is stationary 
freq = np.random.normal(0,0.1,size=N) #oscillator freq
phase0 = np.random.normal(0,1,size=N) #starting phase of oscillator 

#phase plots
fig1,ax1 = plt.subplots(nrows=2,ncols=3)
fig2,ax2 = plt.subplots(nrows=2,ncols=3)
fig1.tight_layout()
fig2.tight_layout()
for i in range(2):
    for j in range(3):
        print(i,j)
        t,phase,R,Theta = evolution(K[i][j],phase0,freq,T=T)
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
            ax2[i][j].set_ylabel("R")
        ax1[i][j].legend(loc="upper center")
        ax1[i][j].set_title("K = {}, N={}".format(K[i][j],N))
        ax1[i][j].set_ylim([-np.pi,np.pi+2]) #to show the legend the right way 
        ax2[i][j].set_title("K = {}, N={}".format(K[i][j],N))
        ax2[i][j].set_ylim([0,1.3])


# Step 2:

#plotting R,Theta for different values of K
K = np.linspace(0,18,100) #take K=0-1
N = 1000 #adapt this to get different amount of neurons 
T = 50 #Take T=100 for Step 2 to see that T=50 most of the time r is stationary 
freq = np.random.normal(0,0.1,size=N) #oscillator freq
phase0 = np.random.normal(0,1,size=N) #starting phase of oscillator 

R_array = [] 
for i in range(len(K)):
    print(i)
    t,phase,R,Theta = evolution(K[i],phase0,freq,T=T)
    R_array.append(np.average(R[-30:-1])) #look at average of last 30 elements of R (because at the end times, R is stationary)
R_array = np.array(R_array) 

fig,ax = plt.subplots(nrows=1,ncols=1)
ax.plot(K,R_array)
ax.set_xlabel('K')
ax.set_ylabel('R')
ax.set_title('stationary R behaviour')


#Step 3:
K = np.linspace(0,15,100) #Step2: take K=0-1
N = 1000 #adapt this to get different amount of neurons 
T = 100 #Take T=70 such that you can average between T=30 to T=70 for the phase difference of a neuron and Theta
freq = np.random.normal(0,0.1,size=N) #oscillator freq
phase0 = np.random.normal(0,1,size=N) #starting phase of oscillator 

phase_array = np.zeros([len(K),N,T*10]) #T*10 is T/h and need to be adapted if h is adapted
Theta_array = np.zeros([len(K),T*10])
for i in range(len(K)):
    print(i)
    t,phase,R,Theta = evolution(K[i],phase0,freq,T=T)
    phase_array[i] = phase
    Theta_array[i] = Theta

difference = np.empty([len(K),N])
for k in range(len(K)):
    Theta = Theta_array[k]
    for n in range(N):
        phase = phase_array[k][n]
        difference[k][n] = np.abs(np.average(phase[-100:-1]%(2*np.pi)-Theta[-100:-1]%(2*np.pi)))
threshold = 0.2 * 2*np.pi
fraction = np.average(np.where(difference<threshold,1,0)[:],axis=1)

fig,ax = plt.subplots(nrows=1,ncols=1)
ax.plot(K,fraction)

#   The amount locked oscillators goes to one for bigger K. This is because if w_i < K*r then 
#   they lock. This means that if K grows the amount of neurons with w_i < K*r grows. 
