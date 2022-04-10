from re import A
import numpy as np 
import matplotlib.pyplot as plt
from tqdm import trange,tqdm
from mpl_toolkits.mplot3d import art3d
from scipy.interpolate import interp2d

def create_circle(n,r,include_border=False): 
    """ Create circle of radius r with approximately n points. Return array with 
    coordinates of these points """
    k = 2*int(np.sqrt(n/np.pi))+2 # this gives the amount of points on a diameter to get approximately n points
    x = np.linspace(-r,r,k)
    y = np.linspace(-r,r,k)
    X,Y = np.meshgrid(x,y)
    if include_border:
        x_in_circle = x[np.where(X**2+Y**2<=r**2)[0]]
        y_in_circle = y[np.where(X**2+Y**2<=r**2)[1]]
    else:
        x_in_circle = x[np.where(X**2+Y**2<r**2)[0]]
        y_in_circle = y[np.where(X**2+Y**2<r**2)[1]]
    points = np.transpose(np.array([x_in_circle,y_in_circle]))
    return points

def show_points(points,r):
    """ Takes points that are supposed to be in circle with radius r and
    visualizes the points inside a cricle. """
    fig,ax = plt.subplots()
    circle = plt.Circle((0,0),r,fill=False)
    ax.scatter(points[:,0],points[:,1],s=0.5)
    ax.add_patch(circle)
    ax.set_xlim([-r-r/10,r+r/10])
    ax.set_ylim([-r-r/10,r+r/10])
    plt.show()

def f(ix,iy,px,py,sigma):
    """ Takes coordinate(s) of place cells and space grid and return the firing rate. The shape it returns depends on the input. """
    if type(ix)==np.ndarray and type(px)==np.ndarray:
        sum1 = np.transpose(np.outer(np.ones(len(px)),np.ones(len(ix)))*ix) - np.outer(np.ones(len(ix)),np.ones(len(px)))*px
        sum2 = np.transpose(np.outer(np.ones(len(py)),np.ones(len(iy)))*iy) - np.outer(np.ones(len(iy)),np.ones(len(py)))*py
        firing_rate = np.exp(-(sum1**2+sum2**2)/(2*sigma**2))
    else:
        firing_rate = np.exp(-((ix-px)**2+(ix-py)**2)/(2*sigma**2))
    return firing_rate

def get_C(place_cells,px,py,w,sigma): 
    """ Return C(p) for all points p given some w. """ 
    firing_rate = f(place_cells[:,0],place_cells[:,1],px,py,sigma)
    if type(px)==np.ndarray:
        C = np.matmul(np.transpose(firing_rate),w)
    else:
        C = np.dot(firing_rate,w)
    return C # otherwise infinities can arise

def get_a(place_cells,px,py,z,sigma): 
    """ Return a(p) for all points p, all directions j given some w. """
    firing_rate = f(place_cells[:,0],place_cells[:,1],px,py,sigma)
    a = np.matmul(z,firing_rate)
    return a # otherwise infinities can arise

def get_R(Rx,Ry,Rr,px,py):
    R = np.where((px-Rx)**2+(py-Ry)**2<Rr**2,1,0)
    return R

def get_neighbours(coord,px,py):
    unique_values = np.unique([px,py])
    x,y = np.where(unique_values==coord[0])[0][0],np.where(unique_values==coord[1])[0][0]
    neighbours = np.zeros([8,2])*np.nan
    for i,(a,b) in enumerate([(0,1),(1,1),(1,0),(1,-1),(0,-1),(-1,-1),(-1,0),(-1,1)]):
        if not(x+a>0 and y+b>0 and x+a<len(unique_values) and y+b<len(unique_values)):
            neighbours[i] = np.nan
        elif unique_values[x+a]**2+unique_values[y+b]**2>r**2:
            neighbours[i] = np.nan
        elif a!=0 or b!=0:
            neighbours[i] = unique_values[x+a],unique_values[y+b]
    return neighbours

def contour_plot(px,py,z,r,Rx,Ry,Rr,D3=False): 
    if D3:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.plot_trisurf(px,py,z,cmap="rainbow")
        contour = ax.tricontourf(px,py,z,cmap="rainbow", offset=0)
        # ax.scatter(place_cells[:,0],place_cells[:,1],0,s=10,color='r') 
        circle = plt.Circle((0,0),r,fill=False,color='r')
        ax.add_patch(circle)
        art3d.pathpatch_2d_to_3d(circle,z=0,zdir='z')
        circle2 = plt.Circle((Rx,Ry),Rr,fill=False,color='blue')
        ax.add_patch(circle2)
        art3d.pathpatch_2d_to_3d(circle2,z=0,zdir='z')
    else:
        fig,ax = plt.subplots()
        contour = ax.tricontourf(px,py,z,cmap="rainbow")
        # ax.scatter(place_cells[:,0],place_cells[:,1],s=10,color='r') 
        circle = plt.Circle((0,0),r,fill=False,color='r')
        ax.add_patch(circle)
        circle2 = plt.Circle((Rx,Ry),Rr,fill=False,color='blue')
        ax.add_patch(circle2)
    fig.colorbar(contour)
    plt.show()

def quiver_plot(px,py,a,Rx,Ry,Rr,stream=False): ######## DEZE FUNCTIE IS NOT NIET AF ############
    fig,ax = plt.subplots()
    theta = np.outer(np.arange(8)*np.pi/4,np.ones(len(place_cells)))
    x,y = get_a(place_cells,place_cells[:,0],place_cells[:,1],z,sigma)*np.array([np.cos(theta),np.sin(theta)])
    ax.quiver(place_cells[:,0],place_cells[:,1],np.sum(x,axis=0),np.sum(y,axis=0),color='r')
    circle = plt.Circle((0,0),r,fill=False,color='r')
    ax.add_patch(circle)
    circle2 = plt.Circle((Rx,Ry),Rr,fill=False,color='blue')
    ax.add_patch(circle2)
    if stream:
        xi = np.linspace(px.min(), px.max(), px.size)
        yi = np.linspace(py.min(), py.max(), py.size)
        uCi = interp2d(px, py, np.sum(x,axis=0))(xi, yi)
        vCi = interp2d(px, py, np.sum(y,axis=0))(xi, yi)
        ax.streamplot(xi, yi, uCi, vCi)
    plt.show()

def path_plot(p_path,trial,space_grid,w,z,t,i,r,Rx,Ry,Rr,sigma):
    fig,ax = plt.subplots() 
    platform_circle = plt.Circle((Rx,Ry),Rr,fill=False,color='r') 
    ax.add_patch(platform_circle) 
    space_circle = plt.Circle((0,0),r,fill=False,color='b') 
    ax.add_patch(space_circle) 
    ax.plot(p_path[t+1:i,0],p_path[t+1:i,1],color='black') 
    contour = ax.tricontourf(space_grid[:,0],space_grid[:,1],get_C(place_cells,space_grid[:,0],space_grid[:,1],w,sigma),cmap="rainbow")
    theta = np.outer(np.arange(8)*np.pi/4,np.ones(len(place_cells)))
    x,y = get_a(place_cells,place_cells[:,0],place_cells[:,1],z,sigma)*np.array([np.cos(theta),np.sin(theta)])
    ax.quiver(place_cells[:,0],place_cells[:,1],np.sum(x,axis=0),np.sum(y,axis=0),color='r')
    ax.set_title('trial {}'.format(trial)) 
    fig.colorbar(contour)
    plt.show() 

def rodent_swim(place_cells,space_grid,w,z,T,gamma,Rx,Ry,Rr,kw,kz,sigma):
    p_start_x = space_grid[np.where(space_grid[:,0]==np.max(space_grid[:,0]))]
    p_start = p_start_x[np.where(np.abs(p_start_x[:,1])==np.min(np.abs(p_start_x[:,1])))][0]
    p = p_start.copy()
    p_path = np.zeros([T,2])
    p_path[0] = p
    trial = 1
    t = 0
    iterations = trange(1,T,desc="trial {}".format(trial))
    count = 0 
    path_plot(p_path,trial,space_grid,w,z,t,1,r,Rx,Ry,Rr,sigma)
    for i in iterations:
        neighbours = get_neighbours(p,space_grid[:,0],space_grid[:,1])
        a = get_a(place_cells,p[0],p[1],z,sigma)*(neighbours[:,0]*0+1)
        C = get_C(place_cells,p[0],p[1],w,sigma)
        R = get_R(Rx,Ry,Rr,p[0],p[1])
        firing_rate = f(place_cells[:,0],place_cells[:,1],p[0],p[1],sigma)
        P_neighbours_nan = np.exp(2*a)/np.nansum(np.exp(2*a))
        P_neighbours = np.where(np.isnan(P_neighbours_nan),0,P_neighbours_nan) # if at boundary this direction cannot be chosen and P=0
        if np.sum(P_neighbours)<0.9:
            print(a,P_neighbours)
        direction = np.random.choice(np.arange(len(P_neighbours)),p=P_neighbours) # pick a direction from the probability distribution
        p = neighbours[direction]
        p_path[i] = p
        if R!=1:
            C_new = get_C(place_cells,p[0],p[1],w,sigma)
        else:
            C_new = 0
        p_path[i] = p
        delta = R+gamma*C_new-C
        dw = delta*firing_rate
        g = np.zeros([8,N])
        g[direction,:] = 1
        dz = delta*firing_rate*g
        if i%1000==0:
            iterations.set_description("trial {}, δ={:.0e}".format(trial,np.abs(delta)))
            iterations.refresh()
        w += kw*dw
        z += kz*dz
        if R==1:
            if trial%20==0:
                path_plot(p_path,trial,space_grid,w,z,t,i,r,Rx,Ry,Rr,sigma)
            p = p_start.copy() 
            t = i
            trial += 1
            count = 0
            iterations.set_description("trial {}, δ={:.0e}".format(trial,np.abs(delta)))
            iterations.refresh()
        count += 1
        if count%10000==0:
            path_plot(p_path,trial,space_grid,w,z,t,i,r,Rx,Ry,Rr,sigma)
            print(P_neighbours)
            if count > 100000:
                print('P',P_neighbours)
                print('Reset')
                p = p_start.copy()
                count = 0
                t = i
    path_plot(p_path,trial,space_grid,w,z,t,t,r,Rx,Ry,Rr,sigma)
    return w,z

r = 2
n = 493 # approximately  points
T = 100000
sigma = 0.2
gamma = 0.9
kw = 0.1 # is learning rate for w and z
kz = 0.1
Rx,Ry,Rr = -0.5,0.5,0.2
place_cells = create_circle(n,r,T)
space_grid = create_circle(n*2,r,include_border=True)
N = len(place_cells) # there are exactly N place_cells
# w0 = np.random.normal(0,1,size=N)
# z0 = np.random.normal(0,1,size=[8,N])
w0 = np.zeros(N)
z0 = np.zeros([8,N])

w,z = rodent_swim(place_cells,space_grid,w0,z0,T,gamma,Rx,Ry,Rr,kw,kz,sigma)

C = get_C(place_cells,space_grid[:,0],space_grid[:,1],w,sigma)
contour_plot(space_grid[:,0],space_grid[:,1],C,r,Rx,Ry,Rr,D3=True)
# a = get_a(place_cells,place_cells[:,0],place_cells[:,1],z,sigma)
# quiver_plot(place_cells[:,0],place_cells[:,1],a,Rx,Ry,Rr,stream=False