import numpy as np
import matplotlib.pyplot as plt

#1a
a0 = 0.5*np.pi
a1 = 1.5*np.pi
b0 = 0*np.pi
b1 = 3*np.pi
L1 = 1
L2 = 1
t = np.linspace(a0,a1,1000)

x_ = L1*np.cos(t)
y_ = L1*np.sin(t) 
x = L1*np.cos(t)+L2*np.cos((b1-b0)/(a1-a0)*t+(b0*a1-a0*b1)/(a1-a0))
y = L1*np.sin(t)+L2*np.sin((b1-b0)/(a1-a0)*t+(b0*a1-a0*b1)/(a1-a0))

plt.plot(x,y,color='b')
plt.plot(x_,y_,color='r')


#1b
x0, y0 = 1, 3
x1, y1 = 4, 1
L1, L2 = 1, 0.5
t = np.linspace(x0,x1,1000)

x = t
y = (y1-y0)/(x1-x0)*t+(y0*x1-x0*y1)/(x1-x0)

plt.plot(x,y,color='b')


#2d/e
#Je kan hier x_0,x_T,y_0,y_T aanpassen zonder dat de vorm verandert. 
x_0, x_T = 0, 10
y_0, y_T = 0, 10
T = 1.0
t = np.linspace(0,int(T),int(T*20))
v = np.sqrt((x_T-x_0)**2+(y_T-y_0)**2)/T*(-30*t**2+60*t**3-30*t**4)
a = np.sqrt((x_T-x_0)**2+(y_T-y_0)**2)/T*(-60*t**2+180*t**3-120*t**4)

plt.plot(t,v,label='v')
plt.plot(t,a,label='a')
plt.title('(x0,y0)=(4,5),    (xT,yT)=(1,6)')
plt.xlabel('t')
plt.ylabel('v')
