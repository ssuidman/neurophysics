import numpy as np
import matplotlib.pyplot as plt


# exercise 1.1
def ex1():
    r0 = 100
    T = 10
    dt = 0.001
    t = np.linspace(0,T,round(T/dt)) 
    spikes = np.zeros(len(t))
    for i in range(len(t)):
        spikes[i] = np.random.poisson(r0*dt)
    ISI = np.diff(np.where(spikes==1))[0]*dt

    CV = np.std(ISI)/np.mean(ISI)
    spikes_reshaped = np.sum(spikes.reshape(round(0.1/dt),round(len(spikes)/round(0.1/dt))),0)
    fano = np.std(spikes_reshaped)/np.mean(spikes_reshaped)

    print("CV:   ",round(CV,3))
    print("fano: ",round(fano,3))

    fig,ax = plt.subplots()
    ax.hist(ISI,bins=30)
    ax.set_xlim(0,0.05)
    ax.set_title("ISI histogram")
    ax.set_xlabel('t(s)')
    fig.show()



def ex2():
    r0 = 100
    r = r0 
    T = 10 
    dt = 0.01 
    N = round(T/dt)
    t = np.linspace(0,T,N) 

    CV_list = []
    ISI_list = []
    tau_values = (np.linspace(1,20,1000))/1000 
    if 0.01 not in tau_values:
        tau_values = np.sort(np.append(tau_values,0.01)) # tau=0.01 must be there
    exp_rate = False
    for tau in tau_values:
        t_exp = 0
        spikes = np.zeros(N)
        for i,j in enumerate(t):
            spike = np.random.poisson(r*dt)
            spikes[i] = spike
            if exp_rate:
                r = r0*(1-np.exp(-t_exp/tau)) 
            if spike==1:
                r = 0
                exp_rate = True
                t_exp = 0
            t_exp += dt
        ISI = np.diff(np.where(spikes==1))[0]*dt 
        ISI_list.append(ISI)
        CV = np.std(ISI)/np.mean(ISI) 
        CV_list.append(CV)
        if tau == 0.01:
            x = spikes
            spikes_reshaped = np.array(np.array_split(spikes,round(T*10)))
            spikes_summed = np.array([np.sum(i) for i in spikes_reshaped])
            fano = np.std(spikes_summed)/np.mean(spikes_summed) 

    print("fano: ",round(fano,3))

    # coefficient of variation as function of τ_ref
    fig,ax = plt.subplots()
    ax.set_title("CV for different τ_ref")
    ax.plot(tau_values*1000,CV_list)
    ax.set_xlabel("τ_ref(ms)")
    ax.set_ylabel("CV")
    fig.show()

    #ISI interval plot for different τ_ref 
    fig,ax = plt.subplots(nrows=2,ncols=2)
    bin_len = 15
    x_min, x_max, y_min, y_max = 0,100 , 0,150

    ax[0][0].hist(ISI_list[0]*1000,bins=bin_len)
    ax[0][0].set_xlabel("t(ms)")
    ax[0][0].set_title("τ_ref={}ms".format(round(tau_values[0]*1000)))
    ax[0][0].set_xlim([x_min,x_max])
    ax[0][0].set_ylim([y_min,y_max])

    ax[0][1].hist(ISI_list[int(0.25*len(tau_values))]*1000,bins=bin_len)
    ax[0][1].set_xlabel("t(ms)")
    ax[0][1].set_title("τ_ref={}ms".format(round(tau_values[int(0.25*len(tau_values))]*1000)))
    ax[0][1].set_xlim([x_min,x_max])
    ax[0][1].set_ylim([y_min,y_max])

    ax[1][0].hist(ISI_list[int(0.5*len(tau_values))]*1000,bins=bin_len)
    ax[1][0].set_xlabel("t(ms)")
    ax[1][0].set_title("τ_ref={}ms".format(round(tau_values[int(0.5*len(tau_values))]*1000)))
    ax[1][0].set_xlim([x_min,x_max])
    ax[1][0].set_ylim([y_min,y_max])

    ax[1][1].hist(ISI_list[int(0.75*len(tau_values))]*1000,bins=bin_len)
    ax[1][1].set_xlabel("t(ms)")
    ax[1][1].set_title("τ_ref={}ms".format(round(tau_values[int(0.75*len(tau_values))]*1000)))
    ax[1][1].set_xlim([x_min,x_max])
    ax[1][1].set_ylim([y_min,y_max])

    fig.suptitle("ISI histograms")
    fig.tight_layout()
    fig.show()



# exercise 1.3
def ex3():
    r0 = 100
    T = 10
    dt = 0.001
    t = np.linspace(0,T,round(T/dt)) 
    N = len(t)
    spikes = np.zeros(N)

    #constant firing rate
    for i,j in enumerate(t):
        spikes[i] = np.random.poisson(r0*dt)
    ISI_const = np.diff(np.where(spikes==1))[0]*dt

    #exponential firing rate
    r = r0
    tau = 0.01
    t_exp = 0
    exp_rate = False

    for i,j in enumerate(t):
        spike = np.random.poisson(r*dt)
        spikes[i] = spike
        if exp_rate:
            r = r0*(1-np.exp(-t_exp/tau)) 
        if spike==1:
            r = 0
            exp_rate = True
            t_exp = 0
        t_exp += dt
    ISI_exp = np.diff(np.where(spikes==1))[0]*dt

    #sinusoidal firing rate
    r = r0*(1+np.cos(2*np.pi*t/0.025))

    for i,j in enumerate(t):
        spike = np.random.poisson(r[i]*dt)
        spikes[i] = spike
    ISI_sin = np.diff(np.where(spikes==1))[0]*dt


    #plots
    fig,ax = plt.subplots(ncols=3)

    ax[0].hist(ISI_const*1000,bins=30)
    ax[0].set_xlim(0,50)
    ax[0].set_title("constant rate")
    ax[0].set_xlabel('t(ms)')

    ax[1].hist(ISI_exp*1000,bins=30)
    ax[1].set_xlim(0,50)
    ax[1].set_title("exponential rate")
    ax[1].set_xlabel('t(ms)')

    ax[2].hist(ISI_sin*1000,bins=30)
    ax[2].set_xlim(0,50)
    ax[2].set_title("sinusoidal rate")
    ax[2].set_xlabel('t(ms)')

    fig.suptitle("ISI histogram")
    fig.tight_layout()
    fig.show()




# exercise 5.3
def ex5():
    
    Ie_array = np.linspace(0,10,100)*10**(-9)
    r_book_array = np.zeros(len(Ie_array))
    r_calc_array = np.zeros(len(Ie_array))
    
    for k,I_e in enumerate(Ie_array):
        #set constants
        E_L = -70*10**(-3)
        R_m = 10*10**6
        tau_m = 10*10**(-3)
        V_th = -54*10**(-3)
        V_reset = -80*10**(-3)
        t_pulse = 300*10**(-3)
        t_simulation = 500*10**(-3)
        dt = 0.0001

        #create loop for dV/dt
        t = np.linspace(0,t_simulation,int(t_simulation/dt))
        V = np.zeros(len(t))
        V[0] = E_L
        V[-1] = E_L
        for i,j in enumerate(t):
            if j<t_pulse:
                dV = (E_L - V[i-1] + R_m*I_e)/tau_m*dt
            else:
                dV = (E_L - V[i-1])/tau_m*dt
            V[i] = V[i-1] + dV
            if V[i] > V_th:
                V[i] = V_reset

        #calculate firing rate r
        r_calc = (np.mean(np.diff(np.where(V==V_reset))*dt))**-1
        r_book = (tau_m*np.log((R_m*I_e+E_L-V_reset)/(R_m*I_e+E_L-V_th)))**-1

        r_calc_array[k] = r_calc
        r_book_array[k] = r_book
    fig,ax = plt.subplots()
    ax.plot(Ie_array*10**9,r_calc_array,label="calculated")
    ax.plot(Ie_array*10**9,r_book_array,label="from book")
    ax.set_xlabel("I(nA)")
    ax.set_ylabel("r(Hz)")
    ax.legend()
    fig.show()
    



ex1()
ex2()
ex3()
ex5()