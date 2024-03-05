import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

## opdr 2a)
def f(x,y):
    return -2*x-y

def second_order(x0,y0,h,f,xlim=[0,1]):
    x = [x0]
    y = [y0]
    xn = x[0]
    yn = y[0]
    while xn>=xlim[0] and xn<xlim[1]:
        yn_tilde = yn + h*f(xn,yn)
        yn += 0.5*h*(f(xn,yn)+f(xn+h,yn_tilde))
        xn += h
        x.append(xn)
        y.append(yn)
    x = np.array(x)
    y = np.array(y)
    return x,y

def fourth_order(x0,y0,h,f,xlim=[0,1]):
    x = [x0]
    y = [y0]
    xn = x[0]
    yn = y[0]
    while xn>=xlim[0] and xn<xlim[1]:
        k1 = h*f(xn,yn)
        k2 = h*f(xn+0.5*h,yn+0.5*k1)
        k3 = h*f(xn+0.5*h,yn+0.5*k2)
        k4 = h*f(xn+h,yn+k3)
        yn += 1/6*(k1+2*k2+2*k3+k4)
        xn += h
        x.append(xn)
        y.append(yn)
    x = np.array(x)
    y = np.array(y)
    return x,y

x2,y2 = second_order(0,-1,0.1,f)
x4,y4 = fourth_order(0,-1,0.1,f)
x = np.linspace(0,1,1001)
y_true = -3*np.exp(-x)-2*x+2

## opdr 2b)
y_odeint = odeint(f,-1,x,tfirst=True).transpose()[0]

plt.plot(x2,y2,label='second order')
plt.plot(x4,y4,label='fourth order')
#plt.plot(x,y_true,'r')
plt.plot(x,y_odeint,label='odeint')
plt.legend()



## opdr 2d) 

def pendulum(x0,y0,h,f,xlim=[0,1],I=0):
    y1 = [y1_0]
    y2 = [y2_0]
    y1_n = y1[0]
    y2_n = y2[0]
    t = 0
    while t<t_end:
        k1 = h*f1()
        k2 = h*f1(xn+0.5*h,yn+0.5*k1)
        k3 = h*f1(xn+0.5*h,yn+0.5*k2)
        k4 = h*f1(xn+h,yn+k3)
        yn += 1/6*(k1+2*k2+2*k3+k4)


    x = [x0]
    y = [y0]
    xn = x[0]
    yn = y[0]
    while xn>=xlim[0] and xn<xlim[1]:
        k1 = h*f(xn,yn)
        k2 = h*f(xn+0.5*h,yn+0.5*k1)
        k3 = h*f(xn+0.5*h,yn+0.5*k2)
        k4 = h*f(xn+h,yn+k3)
        yn += 1/6*(k1+2*k2+2*k3+k4)
        xn += h
        x.append(xn)
        y.append(yn)
    x = np.array(x)
    y = np.array(y)
    return x,y

def y2_prime(y1,y2):
    return (1-np.sin(y1))/y2

#y1 = np.linspace(-np.pi/2-1,np.pi/2+1)
x_pend,y_pend = fourth_order(-np.pi/2,0,0.1,y2_prime,xlim=[-np.pi/2-1,np.pi/2+1])
