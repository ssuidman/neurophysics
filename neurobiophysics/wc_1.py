import numpy as np
import matplotlib.pyplot as plt

def run(N,o):
    x = 0
    y = 0
    path_x = [x]
    path_y = [y]
    for i in range(N):
        x += np.random.normal(0,o)
        y += np.random.normal(0,o)
        path_x.append(x)
        path_y.append(y)
    path_x = np.array(path_x)
    path_y = np.array(path_y)
    #plt.scatter(0,0,color='r',s=20,zorder=1,marker='x')
    #plt.plot(path_x,path_y,linewidth=0.2,zorder=0)
    #plt.show()
    return path_x,path_y

steps = 1000
simulations = 1000
D = []

for sigma in [0.2,0.5,1.0,1.5]:
    length = []
    for sim in range(simulations):
        x_sim,y_sim = run(steps,sigma)
        length_sim = np.sqrt(x_sim**2+y_sim**2)
        length.append(length_sim)
        print(sim,sigma)
    length = np.average(np.array(length),axis=0)
    D_sigma = length**2/(2*1)
    D.append(D_sigma)

D_slope = np.nanmean(D/t,axis=1)
print(D_slope)

t = np.linspace(0,1000,1001)
plt.plot(t,D[0],label='σ=0.2')
plt.plot(t,D[1],label='σ=0.5')
plt.plot(t,D[2],label='σ=1.0')
plt.plot(t,D[3],label='σ=1.5')
plt.legend()
plt.title('Diffusion law')
plt.xlabel('t')
plt.ylabel('D')
plt.show()

