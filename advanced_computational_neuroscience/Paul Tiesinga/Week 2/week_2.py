import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import newton,fsolve
from tqdm import tqdm 
from scipy.signal import find_peaks
from scipy.interpolate import interp1d,interp2d,UnivariateSpline,CubicSpline

a,b,g,u_T = 2,-1,10,0.5 

def f(u_array,a,b,g,u_T): 
    du_dt = -u_array + 1/(1+np.exp(-g*(a*u_array+b-u_T))) 
    return du_dt 

def f_prime(u_array,a,b,g,u_T): 
    d2u_dt2 = -1 + (g*a*np.exp(-g*(a*u_array+b-u_T)))/(1+np.exp(-g*(a*u_array+b-u_T)))**2 
    return d2u_dt2 

def bifurcation(a_array,b,g,u_T): 
    fig,ax = plt.subplots() 
    for i,a in enumerate(tqdm(a_array)): 
        all_solutions = newton(f,np.arange(0,5,0.1),args=(a,b,g,u_T)) 
        round_solutions = np.round(all_solutions,2) 
        solutions = np.unique(round_solutions) 
        if len(solutions)<=3:
            color_array = np.where(f_prime(solutions,a,b,g,u_T)<0,'blue','red')
            ax.scatter(np.ones(len(solutions))*a,solutions,color=color_array,s=0.5) 
    ax.scatter(1.8882174388,0.9493,color='black',s=50)
    ax.set_ylim(-1,2.5)
    ax.set_xlabel("$\\alpha  \\rightarrow$")
    ax.set_ylabel('u  $\\rightarrow$')
    ax.set_title('Bifurcation diagram')
    plt.show() 

a_array = np.linspace(0,8,1000) 
bifurcation(a_array,b,g,u_T) 


####### opdracht 2 ###########

def F(u,o): 
    f = u**2/(u**2+o**2) 
    return f 

def G(v,o): 
    g = o*np.sqrt(v/(1-v)) 
    return g 

def fixed_points(w11,w22,w12,w21,o,test=False):
    x = np.linspace(0,1,10000)
    y1 = (G(x,o)-w11*x)/w12 
    y2_test = np.linspace(0,1,10000) 
    x2_test = (G(y2_test,o)-w22*y2_test)/w21 
    if find_peaks(np.abs(x2_test))[0]!=[]:
        i = find_peaks(np.abs(x2_test))[0][0]
        j = find_peaks(-np.abs(x2_test))[0][0]

        x2_test1 = x2_test[:i] # segment 1
        y2_test1 = y2_test[:i] # segment 1
        x1 = x[np.argmin(np.abs(x-np.min(x2_test1))):np.argmin(np.abs(x-np.max(x2_test1)))]
        y1_1 = y1[np.argmin(np.abs(x-np.min(x2_test1))):np.argmin(np.abs(x-np.max(x2_test1)))]
        y2_1 = interp1d(x2_test1,y2_test1,bounds_error=False)(x1)

        x2_test2 = x2_test[i:j] # segment 2
        y2_test2 = y2_test[i:j] # segment 2
        x2 = x[np.argmin(np.abs(x-np.min(x2_test2))):np.argmin(np.abs(x-np.max(x2_test2)))]
        y1_2 = y1[np.argmin(np.abs(x-np.min(x2_test2))):np.argmin(np.abs(x-np.max(x2_test2)))]
        y2_2 = interp1d(x2_test2,y2_test2,bounds_error=False)(x2)

        x2_test3 = x2_test[j:-1] # segment 3 
        y2_test3 = y2_test[j:-1] # segment 3 
        x3 = x[np.argmin(np.abs(x-np.min(x2_test3))):np.argmin(np.abs(x-np.max(x2_test3)))] 
        y1_3 = y1[np.argmin(np.abs(x-np.min(x2_test3))):np.argmin(np.abs(x-np.max(x2_test3)))] 
        y2_3 = interp1d(x2_test3,y2_test3,bounds_error=False)(x3) 

        x0,y0 = [],[] 
        x0.append(x1[find_peaks(-np.abs(y2_1-y1_1))[0]]) 
        y0.append(y1_1[find_peaks(-np.abs(y2_1-y1_1))[0]]) 
        x0.append(x2[find_peaks(-np.abs(y2_2-y1_2))[0]]) 
        y0.append(y1_2[find_peaks(-np.abs(y2_2-y1_2))[0]]) 
        x0.append(x3[find_peaks(-np.abs(y2_3-y1_3))[0]]) 
        y0.append(y1_3[find_peaks(-np.abs(y2_3-y1_3))[0]]) 
        x0 = np.concatenate(x0) 
        y0 = np.concatenate(y0) 
    else: 
        y2 = interp1d(x2_test,y2_test)(x) 
        x0 = x[find_peaks(-np.abs(y1-y2))[0]]
        y0 = y1[find_peaks(-np.abs(y1-y2))[0]]
        # x0,y0 = 0,0
    if test:
        return x0,y0,x,y1,x2_test,y2_test
    else:
        return x0,y0

def find_bifurcation(w11,w22,w12,w21):
    o_array = np.linspace(0.1,0.4,1000) 
    fp_length = np.zeros(len(o_array))
    for i,o in enumerate(tqdm(o_array)):
        x0,y0 = fixed_points(w11,w22,w12,w21,o)
        fp_length[i] = len(x0)
    fp_length = np.where(fp_length>3,np.nan,fp_length)
    return fp_length,o_array

w11,w22,w12,w21,o = 0.3,0.2,0.2,0.3,0.2 
x0,y0,x,y1,x2_test,y2_test = fixed_points(w11,w22,w12,w21,o,test=True) 

fig,ax = plt.subplots() 
ax.plot(x,y1) 
ax.plot(x2_test,y2_test) 
ax.scatter(x0,y0) 
ax.set_xlim([-0.01,1]) 
ax.set_ylim([-0.4,1]) 
ax.set_title('$\sigma$=0.2')
ax.set_xlabel('u1')
ax.set_xlabel('u2')
plt.show()

fp_length,o_array = find_bifurcation(w11,w22,w12,w21) 
bifurcations = o_array[np.where(np.diff(fp_length)!=0)] 
print("The changes in the length of fixed points are at:",bifurcations)


