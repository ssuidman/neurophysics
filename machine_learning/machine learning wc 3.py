import numpy as np
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import diagonal

#set some variables
N=10
e = 0.01
P_training = 5
P_test = 10000
x_training = np.random.choice([-1,1],size=(P_training,N))
x_test = np.random.choice([-1,1],size=(P_test,N))

#create teacher 
w_teacher = np.random.normal(0,1,size=N)
y_teacher = np.sign(np.matmul(x_training,w_teacher))

#learning
indices = np.dot(w,x[1])
w*x[i]*y > 0 dan niks
w*x[i]*y < 0 dan w = w + e*x[i]*y

