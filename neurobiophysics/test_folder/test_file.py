import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.get_backend()

points = np.array([[0,0],[1,0],[1/2,np.sqrt(3/4)]])

def f(x,y,N,points):
    x_points = points[:,0] 
    y_points = points[:,1] 
    x_list = [x] 
    y_list = [y] 
    for i in range(N): 
        if i in [int(a*N/100) for a in range(1,101)]: 
            print(str(int(i/N*100))+'%') 
        j = np.random.randint(len(points)) 
        x = x_points[j]*0.5 + x*0.5 
        y = y_points[j]*0.5 + y*0.5 
        x_list.append(x) 
        y_list.append(y) 
    return(x_list,y_list)      

    x_lijst,y_lijst = f(0,0.5,1000000,points) 
    plt.scatter(x_lijst,y_lijst,s=0.5) 
    plt.show() 
