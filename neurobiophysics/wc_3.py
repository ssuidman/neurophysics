import numpy as np
import matplotlib.pyplot as plt

#1a
phi0, phi1 = 0*np.pi, 1.5*np.pi
theta0, theta1 = 0*np.pi, 1*np.pi
dphi = phi1-phi0
dtheta = theta1-theta0
L1, L2 = 1, 0.5 
t = np.linspace(-10,10,1000)
u = 1/(1+np.exp(-t))

x = L1*np.cos(phi0+dphi*u) + L2*np.cos((phi0+theta0)+(dphi+dtheta)*u)
y = L1*np.sin(phi0+dphi*u) + L2*np.sin((phi0+theta0)+(dphi+dtheta)*u)
x_ = L1*np.cos(phi0+dphi*u) 
y_ = L1*np.sin(phi0+dphi*u)


plt.plot(x,y,color='b',label='hand trajectory')
plt.plot(x_,y_,color='r',label='elbow trajectory')
plt.title('phi∈[0,1.5π] theta∈[0,π]')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()


#1b
phi0, phi1 = 0*np.pi, 1.5*np.pi
theta0, theta1 = 0*np.pi, 1*np.pi
dphi = phi1-phi0
dtheta = theta1-theta0
L1, L2 = 1, 1 
t = np.linspace(-100,100,10000)
dt = t[1]-t[0]
u = 1/(1+np.exp(-t))

vx = -L1*dphi*np.sin(phi0+dphi*u)*u*(1-u) - L2*(dphi+dtheta)*np.sin((phi0+theta0)+(dphi+dtheta)*u)*u*(1-u)
vy = L1*dphi*np.cos(phi0+dphi*u)*u*(1-u) + L2*(dphi+dtheta)*np.cos((phi0+theta0)+(dphi+dtheta)*u)*u*(1-u)
vphi = dphi*u*(1-u)
vtheta = dtheta*u*(1-u)

fig, ax = plt.subplots(nrows=2,ncols=2)
ax[0,0].plot(t,vx,label='x',color='b')
ax[1,0].plot(t,vy,label='y',color='r')
ax[0,1].plot(t,vphi,label='phi',color='g')
ax[1,1].plot(t,vtheta,label='theta',color='k')
ax[1,0].set_xlabel('t')
ax[1,1].set_xlabel('t')
ax[0,0].set_ylabel('v')
ax[1,0].set_ylabel('v')
ax[0,0].set_xlim(-10,10)
ax[1,0].set_xlim(-10,10)
ax[0,1].set_xlim(-10,10)
ax[1,1].set_xlim(-10,10)
ax[0,0].legend()
ax[1,0].legend()
ax[0,1].legend()
ax[1,1].legend()
fig.suptitle('velocity profiles')


#2d/e
#Je kan hier x_0,x_T,y_0,y_T aanpassen zonder dat de vorm verandert. 
x_0, x_T = 0, 10
y_0, y_T = 0, 10
T = 1.0
t = np.linspace(0,int(T),int(T*100))
v = np.sqrt((x_T-x_0)**2+(y_T-y_0)**2)/T*(-30*t**2+60*t**3-30*t**4)
a = np.sqrt((x_T-x_0)**2+(y_T-y_0)**2)/T*(-60*t**2+180*t**3-120*t**4)

plt.plot(t,v,label='v')
plt.plot(t,a,label='a')
plt.title('(x0,y0)=(4,5),    (xT,yT)=(1,6)')
plt.xlabel('t')
plt.ylabel('v')
