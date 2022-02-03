import numpy as np
import matplotlib.pyplot as plt
import numpy as np

f = 100
T = 10
k = 100 #for defining how much space should be between two poisson spikes, per second there are f*k timesteps
dt = 1/(f*k) # take time steps k times smaller than histogram
t = np.linspace(0,T,round(T/dt))
spikes = np.zeros(round(T/dt),dtype=int)

for i in range(len(spikes)):
    spikes[i] = np.random.poisson(1/k) # np.sum(spikes)=T*f should be in the end 
ISI = np.diff(np.where(spikes==1))[0]*dt
CV = np.std(ISI)/np.mean(ISI)
spikes_reshaped = np.sum(spikes.reshape(round(0.1/dt),round(len(spikes)/round(0.1/dt))),0)
fano = np.std(spikes_reshaped)/np.mean(spikes_reshaped)
print("CV:   ",round(CV,3))
print("fano: ",round(fano,3))


plt.hist(ISI,bins=30)
plt.xlim(0,0.05)
