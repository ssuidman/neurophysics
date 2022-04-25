import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

T = 0.3 # T = 50 ms 
n = 1000 # is the time grid 
trials = 100 # amount of experiments to be done
p_array = np.linspace(0,0.5,10)
r_sc = np.zeros([len(p_array),trials])

for i,p in enumerate(tqdm(p_array)):
    for j in range(trials):
        t = np.linspace(0,T,n) 
        freq_A = 20*np.abs(np.sin(2*t)) # f_max = 30 Hz 
        freq_B = 30*np.abs(np.sin(3*t)) # f_average = 20 Hz 

        spike_train_A = np.random.poisson(freq_A)
        spike_train_B = np.random.poisson(freq_B)

        spikes_A = np.where(spike_train_A!=0)[0]
        random_spikes_A = np.random.choice(spikes_A,size=int(len(spikes_A)*p),replace=False)
        spike_train_A[random_spikes_A] -= 1
        spike_train_B[random_spikes_A] += 1

        spikes_B = np.where(spike_train_B!=0)[0]
        random_spikes_B = np.random.choice(spikes_B,size=int(len(spikes_B)*p),replace=False)
        spike_train_B[random_spikes_B] -= 1
        spike_train_A[random_spikes_B] += 1

        mean_rate_A = np.mean(spike_train_A)
        mean_rate_B = np.mean(spike_train_B)

        COV_AB = np.sum((spike_train_A-mean_rate_A)*(spike_train_B-mean_rate_B))
        sigma_A = np.sqrt(np.sum((spike_train_A-mean_rate_A)**2))
        sigma_B = np.sqrt(np.sum((spike_train_B-mean_rate_B)**2))
        r_sc_trial = COV_AB/(sigma_A*sigma_B)
        r_sc[i,j] = r_sc_trial


########### GA HIER VERDER, DE GRAFIEKEN KLOPPEN NIET MET 2B UIT DE PAPER!!!!!!!!!! #############
for i in range(10):
    plt.figure()
    plt.hist(r_sc[i])
plt.show()
