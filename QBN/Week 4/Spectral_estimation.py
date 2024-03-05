import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import lfilter,butter,coherence,hann
from tqdm import trange

###########################################################################
################################  1.1.3  ##################################
###########################################################################

a2,f,fs = -0.9,60,1000
w_max = 2*np.pi*f/fs
a1 = 4*a2/(a2-1)*np.cos(w_max)
b,a = np.array([1,1,1]),np.array([1,-a1,-a2])
N = 100000 # must be divisible by 40 

# First create normal GWN from a random signal x
x = np.random.normal(0, 1,N)
Sxx_GWN_0 = np.abs(np.fft.fft(x)/np.sqrt(N))**2
freq = np.fft.fftfreq(N)
freq_GWN = freq[np.where(freq>0)]
Sxx_GWN = Sxx_GWN_0[np.where(freq>0)]

# Now ceate a filtered signal
x_AR2 = lfilter(b,a,x) 
Sxx_AR2_0 = np.abs(np.fft.fft(x_AR2)/np.sqrt(N))**2 
freq_AR2 = freq_GWN.copy()
Sxx_AR2 = Sxx_AR2_0[np.where(freq>0)] 

# Now use Bartlett averaging
x_reshaped = x.reshape(40,int(N/40)) # each x has length 250
x_Bart = lfilter(b,a,x_reshaped,axis=1) 
Sxx_Bart_0 = np.abs(np.fft.fft(x_Bart)/np.sqrt(int(N/40)))**2 
freq_Bart_0 = np.fft.fftfreq(int(N/40)) 
Sxx_Bart_all = Sxx_Bart_0[:,np.where(freq_Bart_0>0)[0]] 
Sxx_Bart = np.mean(Sxx_Bart_all,axis=0)
freq_Bart = freq_Bart_0[np.where(freq_Bart_0>0)[0]] 

# Comparing with Syy
freq_Syy = freq_GWN.copy()
Syy = 1/np.abs(1-a1*np.exp(-1j*2*np.pi*freq_Syy)-a2*np.exp(-2j*2*np.pi*freq_Syy))**2

fig,ax = plt.subplots(nrows=2,ncols=2)

ax[0,0].plot(freq_GWN,Sxx_GWN)
ax[0,0].set_xlabel('f(Hz)')
ax[0,0].set_ylabel('power')
ax[0,0].set_title('GWN')

ax[0,1].plot(freq_AR2,Sxx_AR2,label='AR(2)')
ax[0,1].plot(freq_Bart,Sxx_Bart,label='Bartlett')
ax[0,1].set_xlabel('f(Hz)')
ax[0,1].set_ylabel('power')
ax[0,1].set_title('AR(2)')
ax[0,1].legend()

ax[1,0].plot(freq_Bart,Sxx_Bart*hann(len(Sxx_Bart)))
ax[1,0].set_xlabel('f(Hz)')
ax[1,0].set_ylabel('power')
ax[1,0].set_title('Bartlett with Hann taper')

ax[1,1].plot(freq_Syy,Syy)
ax[1,1].set_xlabel('f(Hz)')
ax[1,1].set_ylabel('power')
ax[1,1].set_title('Syy')

fig.suptitle('Fourier spectra for N={}'.format(N))
fig.tight_layout()
plt.show()

print('The power peak is a little biased towards lower frequencies with Bartlett avering, but barely visible.')


###########################################################################
################################  1.2.1  ##################################
###########################################################################

a2,f,fs = -0.9,60,1000
w_max = 2*np.pi*f/fs
a1 = 4*a2/(a2-1)*np.cos(w_max)
b,a = np.array([1,1,1]),np.array([1,-a1,-a2])
N = 1000 
trials = 1000 

xy_coherence = np.zeros([trials,129]) 
xy_coherence_estim = np.zeros([trials,499])
phase_x = np.zeros([trials,499])
phase_y = np.zeros([trials,499])

for i in trange(trials): 
    x = np.random.normal(0,1,N) 
    x_AR2 = lfilter(b,a,x) 
    noise = np.random.normal(0,1*np.std(x_AR2),N) 

    Sxx_0 = np.abs(np.fft.fft(x_AR2)/np.sqrt(N))**2
    freq = np.fft.fftfreq(N)
    Sxx = Sxx_0[np.where(freq>0)]

    S_noise_0 = np.abs(np.fft.fft(noise)/np.sqrt(N))**2
    S_noise = S_noise_0[np.where(freq>0)]
    
    y_AR2 = x_AR2 + noise
    xy_coherence_trial = coherence(x_AR2,y_AR2) 
    xy_coherence[i] = xy_coherence_trial[1] 
    xy_coherence_estim[i] = Sxx/(S_noise+Sxx)
    fft_x = np.fft.fft(x_AR2)
    fft_x = fft_x[np.where(freq>0)[0]]
    fft_y = np.fft.fft(y_AR2)
    fft_y = fft_y[np.where(freq>0)[0]]

    phase_x[i] = np.angle(fft_x)
    phase_y[i] = np.angle(fft_y)
    # np.abs(np.sum(np.exp(1j*np.angle(fft_x))))/N


freq_coherence = xy_coherence_trial[0] 
freq_coherence_estim = freq[np.where(freq>0)]

fig,ax = plt.subplots(nrows=2) 

ax[0].plot(freq_AR2,Sxx_AR2) 
ax[0].set_ylabel('power') 
ax[0].set_title('Fourier spectrum AR(2)') 

ax[1].plot(freq_coherence,np.mean(xy_coherence,axis=0),label='actual coherence') 
ax[1].plot(freq_coherence_estim,np.mean(xy_coherence_estim,axis=0),label='theoretical coherence') 
ax[1].set_xlabel('f(Hz)') 
ax[1].set_ylabel('power') 
ax[1].set_title('x(t),y(t)-coherence') 

fig.tight_layout() 
plt.show() 

fig,ax = plt.subplots()
ax.plot(np.linspace(0,0.5,int(N/2)-1),np.mean(phase_x,axis=0)-np.mean(phase_y,axis=0))
ax.set_xlabel('f(Hz)')
ax.set_ylabel('phase difference')
plt.show()

PLV = np.abs(np.sum(np.exp(1j*(np.mean(phase_x,axis=0)-np.mean(phase_y,axis=0)))))/len(np.mean(phase_x,axis=0))
print("Phase locking value:",PLV)

