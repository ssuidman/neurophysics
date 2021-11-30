import numpy as np
import matplotlib.pyplot as plt

###############################################################
#   1a  #
###############################################################

def modified_euler(f,v0,u0,h,t0,t1,I):
    v_list = [v0]
    u_list = [u0]
    v = v_list[0]
    u = u_list[0]
    time_range = np.linspace(t0,t1,int((t1-t0)/h))
    for t in time_range:
        v_tilde = v + h*f(v,u,I)[0]
        u_tilde = u + h*f(v,u,I)[1]
        v += 0.5*h*(f(v,u,I)[0]+f(v_tilde,u_tilde,I)[0])
        u += 0.5*h*(f(v,u,I)[1]+f(v_tilde,u_tilde,I)[1])
        if v>35:
            v = -50
            u += 100
        v_list.append(v)
        u_list.append(u)
    return np.array(v_list),np.array(u_list)


def phase_plane(system,x0,y0,a,b=0,c=0,t_range=[0,10],xlim=[-10,10],ylim=[-10,10],variables=["x","y"]): 
    t = np.linspace(t_range[0],t_range[1],(t_range[1]-t_range[0])*100) #set time range
    x,y = modified_euler(f,x0,y0,0.01,t_range[0],t_range[1],a) #creates arrays for x,y
    
    x_range = np.linspace(xlim[0], xlim[1], int((xlim[1]-xlim[0])*6)) #don't use a too big number density (third setting), otherwise the vector field is not clear
    y_range = np.linspace(ylim[0], ylim[1], int((ylim[1]-ylim[0])*6))
    X,Y = np.meshgrid(x_range, y_range)
    dX,dY = system([X,Y],t,a,b,c)

    fig,ax = plt.subplots(nrows=1,ncols=1)
    ax.set_title("phase plane")
    ax.set_xlabel(variables[0])
    ax.set_ylabel(variables[1])   
    ax.set_xlim(xlim[0],xlim[1])
    ax.set_ylim(ylim[0],ylim[1])

    ax.plot(x,y) #plots the trajectory
    ax.streamplot(X,Y,dX,dY, density=1) #plots the streamplot
    #ax.quiver(X,Y,dX,dY,scale=50,scale_units='width',color='blue') #another way to plot the vectorfield
    return fig



#Voor f,system moet je allebei hetzelfde invoeren:
def f(v,u,I):
    dv = 1/100*(0.7*(v+60)*(v+40)-u+I)
    du = 0.03*(-2*(v+60)-u)
    return dv, du

def system(x_y,t,a,b,c):
    x,y = x_y
    dx = 1/100*(0.7*(x+60)*(x+40)-y+a) + 0*x #the 0*x term is such that the shape is always good
    dy = 0.03*(-2*(x+60)-y) + 0*x
    return [dx,dy]

v0 = -50
u0 = 0
tlim = [0,1000]
vlim = [-70,35]
ulim = [-100,300]
h = 0.01
I = 60
v,u = modified_euler(f,v0,u0,h,tlim[0],tlim[1],I)

fig_phase_plane = phase_plane(system,x0=v0,y0=u0,a=I,t_range=tlim,xlim=vlim,ylim=ulim,variables=["v","u"])
plt.figure()
plt.plot(np.linspace(tlim[0],tlim[1],len(v)),v)
plt.title('voltage')
plt.xlabel('time')
plt.ylabel('v')
plt.show()



###############################################################
#   1b  #
###############################################################
"""
-   For I=51 the neurons starts firing if the starting value is V=-50. 

-   A bifurcation happens here. First it is a stable fixed point around 
    V=-52 and no spiking occurs. However at the bifurcation, the system 
    changes to a stable limit cycle where it keeps firing. 

-   The firing rate increases when I is increased. There are for example
    only 3 spikes in 1000 sec for I=55, while this is 13 spikes for I=100. 

"""


###############################################################
#   1c  #
###############################################################


def coupled_neurons(f,N,v0,u0,h,t0,t1,I0,k):
    N = 5
    connection_matrix = [ [0,k,0,0,k] , [k,0,k,0,0] , [0,k,0,k,0] , [0,0,k,0,k] , [k,0,0,k,0] ]
    time_range = np.linspace(t0,t1,int((t1-t0)/h))
    v_list = [[v0 for j in range(len(time_range)+1)] for i in range(N)]
    v_list[0][0] = v0
    u_list = [[0 for j in range(len(time_range)+1)] for i in range(N)]
    v_list[0][0] = u0
    I_list = [[0 for j in range(len(time_range)+1)] for i in range(N)]
    I_list[0][0] = I0
    
    for t in range(len(time_range)):
        for n in range(N):    
            v = v_list[n][t]
            u = u_list[n][t]
            I = I_list[n][t]
            v_tilde = v + h*f(v,u,I)[0]
            u_tilde = u + h*f(v,u,I)[1]
            v += 0.5*h*(f(v,u,I)[0]+f(v_tilde,u_tilde,I)[0]) 
            u += 0.5*h*(f(v,u,I)[1]+f(v_tilde,u_tilde,I)[1]) 
            if v<=35:
                v_list[n][t+1] = v
                u_list[n][t+1] = u
                I_list[n][t+1] += I
            elif v>35:
                v_list[n][t+1] = -50
                u_list[n][t+1] = u+100
                for m in range(N):
                    t_pulse = 0.1 #sec is the duration of a pulse
                    for i in range(int(t_pulse/h)):
                        I_list[m][t+1+i] += connection_matrix[n][m] * I0
        #if I_list[0][t+1] == 0: #does the fist neuron need to get a continuous spike?
        #    I_list[0][t+1] = I0
    return v_list,u_list,I_list #v_list[n][t] gives the voltage of neuron n on time t. This means v_list[n] gives its voltage evolving over time 


v,u,I = coupled_neurons(f,3,-50,0,0.1,0,300,200,0.5)
for i in range(5):
    plt.plot(v[i], label=i)
plt.legend()

plt.eventplot(v)
