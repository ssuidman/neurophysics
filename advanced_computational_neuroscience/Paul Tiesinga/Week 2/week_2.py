import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import newton 
from tqdm import tqdm 

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

w11,w22,w12,w21,o = 0.3,0.2,0.2,0.3,0.2

u1_null1 = np.linspace(0,1,100)
u2_null1 = (G(u1_null1,o)-w11*u1_null1)/w12

u2_null2 = np.linspace(0,1,100)
u1_null2 = (G(u2_null2,o)-w22*u2_null2)/w21

fig,ax = plt.subplots()
ax.plot(u1_null1,u2_null1,label='1')
ax.plot(u1_null2,u2_null2,label='2')
ax.plot(u1_null1,u2_null1-u2_null2,label='difference')
ax.set_xlim([0,2])
ax.set_ylim([0,2])
ax.set_xlabel('$u_1$')
ax.set_ylabel('$u_2$')
ax.set_title('Phase plane')
ax.legend()


print(u1_null1[np.where(np.abs(u2_null1-u2_null2)<10**-2)])
print(u2_null1[np.where(np.abs(u2_null1-u2_null2)<10**-2)])

