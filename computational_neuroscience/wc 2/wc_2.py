import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def vectors(X,Y,t,a,b,c,cart_pol): #In this function you can tell if you're in X,Y or R,Theta. 
    if cart_pol == 'cart':
        dX,dY = system([X,Y],t,a,b,c)
    elif cart_pol == 'pol':
        R = np.sqrt(X ** 2 + Y ** 2)
        Theta = np.arctan(X / Y)
        dR,dTheta = system([R,Theta],t,a,b,c)
        dX = X/np.sqrt(X**2+Y**2)*dR - Y*dTheta
        dY = Y/np.sqrt(X**2+Y**2)*dR + X*dTheta
    return dX,dY



def phase_plane(system,x0,y0,a,b=0,c=0,t_range=[0,10],xlim=[-10,10],ylim=[-10,10],variables=["x","y"]): 
    t = np.linspace(t_range[0],t_range[1],(t_range[1]-t_range[0])*100) #set time range
    x,y = odeint(system,[x0,y0],t,args=(a,b,c)).transpose() #creates arrays for x,y
    
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



def cartesian_plot(system,x0,y0,a,b=0,c=0,t_range=[0,10],xlim=[-10,10],ylim=[-10,10],cart_pol="cart"): 
    t = np.linspace(t_range[0],t_range[1],(t_range[1]-t_range[0])*100) #set time range
    x,y = odeint(system,[x0,y0],t,args=(a,b,c)).transpose() #creates arrays for x,y
    
    x_range = np.linspace(xlim[0], xlim[1], int((xlim[1]-xlim[0])*6)) #don't use a too big number density (third setting), otherwise the vector field is not clear
    y_range = np.linspace(ylim[0], ylim[1], int((ylim[1]-ylim[0])*6))
    X,Y = np.meshgrid(x_range, y_range)
    dX,dY = vectors(X,Y,t,a,b,c,cart_pol)

    fig,ax = plt.subplots(nrows=1,ncols=1)
    ax.set_title("cartesian plot")
    ax.set_xlabel('x')
    ax.set_ylabel('y')   
    ax.set_xlim(xlim[0],xlim[1])
    ax.set_ylim(ylim[0],ylim[1])

    if cart_pol == "cart":
        ax.plot(x,y) #plots the trajectory
    elif cart_pol == "pol":
        ax.plot(x*np.cos(y),x*np.sin(y))
    ax.streamplot(X,Y,dX,dY, density=1) #plots the streamplot
    #ax.quiver(X,Y,dX,dY,scale=50,scale_units='width',color='blue') #another way to plot the vectorfield
    return fig

def system(x_y,t,a,b,c):
    x,y = x_y
    dx = x*(a-x**2) + 0*x #the 0*x term is such that the shape is always good
    dy = 2*np.pi + 0*x
    return [dx,dy]

fig_phase_plane = phase_plane(system,x0=-7.5,y0=7.5,a=2,variables=["r","theta"])
fig_cartesian_plot = cartesian_plot(system,x0=-7.5,y0=7.5,a=2,cart_pol="pol")
fig_phase_plane.show()
fig_cartesian_plot.show()

