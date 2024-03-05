import numpy as np
from numpy import array
from numpy import sum 
from numpy import matmul
from numpy import multiply
from numpy import dot
from numpy import exp
from numpy import log
from numpy import abs
from numpy import min
from numpy import transpose
from numpy.core.function_base import linspace
from numpy.linalg import inv
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy import optimize
import time

#Loading in the MNIST data
file = loadmat('/Users/samsuidman/Desktop/neurophysics/machine_learning/mnistAll.mat') 

#Preparing the data
train_images = file['mnist']['train_images'][0][0]
train_labels = file['mnist']['train_labels'][0][0].transpose()[0]
test_images = file['mnist']['test_images'][0][0]
test_labels = file['mnist']['test_labels'][0][0].transpose()[0]
indices_3 = np.where(train_labels==3)[0]
indices_7 = np.where(train_labels==7)[0]
indices_3_test = np.where(test_labels==3)[0]
indices_7_test = np.where(test_labels==7)[0]
X3 = train_images[:,:,indices_3]
X7 = train_images[:,:,indices_7]
X3_test = test_images[:,:,indices_3_test]
X7_test = test_images[:,:,indices_7_test]
n3 = np.size(X3,2)
n7 = np.size(X7,2)
n3_test = np.size(X3_test,2)
n7_test = np.size(X7_test,2)
X3 = np.reshape(X3,[784,n3])
X7 = np.reshape(X7,[784,n7])
X3_test = np.reshape(X3_test,[784,n3_test])
X7_test = np.reshape(X7_test,[784,n7_test])
X3 = X3/np.max((np.max(np.concatenate(X3)),np.max(np.concatenate(X7))))
X7 = X7/np.max((np.max(np.concatenate(X3)),np.max(np.concatenate(X7))))
X3_test = X3_test/np.max((np.max(np.concatenate(X3_test)),np.max(np.concatenate(X7_test))))
X7_test = X7_test/np.max((np.max(np.concatenate(X3_test)),np.max(np.concatenate(X7_test))))
X3 = (np.insert(X3,0,1,axis=0)).transpose() #shape = (6131,785)
X7 = (np.insert(X7,0,1,axis=0)).transpose() #shape = (6265,785)
X3_test = (np.insert(X3_test,0,1,axis=0)).transpose() #shape = (1010,785)
X7_test = (np.insert(X7_test,0,1,axis=0)).transpose() #shape = (1028,785)
t3 = np.zeros([1,n3])[0] #length = 6131
t7 = np.ones([1,n7])[0] #length = 6265
t3_test = np.zeros([1,n3_test])[0] #length = 1010
t7_test = np.ones([1,n7_test])[0] #length = 1028

#These are the important matrices and arrays in the end. X3 consists of P=6131 images of the number 3. X[3] for example is the 4th image. Each image has 28x28=784 pixels. To make sure that an image doesn't only consists of zeros, before each image is a 1 added. 
X = np.concatenate([X3,X7]) #shape = (12396,785)
t = np.concatenate([t3,t7]) #shape = (12396,)
X_test = np.concatenate([X3_test,X7_test]) #shape = (2038,785)
t_test = np.concatenate([t3_test,t7_test]) #shape = (2038,)





#these are the functions that are defined in the exercise 

def y(w,X): #this function takes a random w and a matrix with N patterns and d dimensions (in our case 785 pixels in total) and returns an array that contains for each pattern the probability that the number is a 3 (or 7)
    y = 1/(1+exp(-dot(X,w)))
    return y 

def E(w,X,t): #Return the error for w,X and the actual value 3,7 (so actual 1,0)
    y0 = y(w,X)
    Ew = -sum(t*log(y0)+(1-t)*log(1-y0))/len(X)
    return Ew

def dE(w,X,t): #gradient for w,X,t
    y0 = y(w,X)
    dE = matmul((y0-t),X)/len(X)
    return dE

def H(w,X): #Hessian for w,X
    y0 = y(w,X) 
    X_y = transpose(multiply(transpose(X),y0*(1-y0))) #multiply each pattern in X by the y value of that pattern. This results in a matrix. 
    H = matmul(transpose(X_y),X)/len(X) #multiply this X*y times X, and sum over all patterns. The new shape is then (785,785)
    return H

def dE_weight_decay(w,X,t,k): #Gradient for the weight decay excercises
    y0 = y(w,X)
    dE = matmul((y0-t),X)/len(X)+k/len(w)*w
    return dE

def H_weight_decay(w,X,k): #Hessian for the weight decay exercises
    y0 = y(w,X) 
    X_y = transpose(multiply(transpose(X),y0*(1-y0))) #multiply each pattern in X by the y value of that pattern. This results in a matrix. 
    H = matmul(transpose(X_y),X)/len(X) + k/len(w)*np.eye(len(w)) #multiply this X*y times X, and sum over all patterns. The new shape is then (785,785)
    return H

X_line_search = X.copy()
t_line_search = t.copy()
def E_line_search(w): #Return the error for w,X and the actual value 3,7 (so actual 1,0)
    y0 = 1/(1+np.exp(-np.dot(X_line_search,w)))
    Ew = -sum(t_line_search*log(y0)+(1-t_line_search)*log(1-y0))/len(X_line_search)
    return Ew

def dE_line_search(w): #gradient for w,X,t
    y0 = 1/(1+np.exp(-np.dot(X_line_search,w)))
    dE = matmul((y0-t_line_search),X_line_search)/len(X_line_search)
    return dE



#These are other helpful functions

def wrong_patterns(w,X,t): #Takes a model w, a dataset X and labels t and return the fraction of misclassified patterns
    y0 = y(w,X)
    classifications = np.where(y0>0.5,1,0)
    wrong_predictions = np.sum(classifications!=t)
    wrong_fraction = wrong_predictions/len(X)
    return wrong_fraction


def visualize(w,X,i): #Takes a model w, a dataset X and a image number and visualizes the image with what the model thinks it represents
    digit = int(np.where(y(w,X[i])>0.5,7,3))
    plt.imshow(X[i][1:].reshape(28,28))
    return digit

def show_results(results,t):
    print('training set:  ' + 'E=' + str(round(results[0][0],3)) + '   wrong patterns: ' + str(round(100*results[0][1],1)) + '%')
    print('test set:      ' + 'E=' + str(round(results[1][0],3)) + '   wrong patterns: ' + str(round(100*results[1][1],1)) + '%')
    print('time passed during learning: ' + str(round(t,1)))

def E_figures(E_training,E_test,title):
    fig,ax = plt.subplots(nrows=1,ncols=1)
    ax.plot(np.log(E_training),label='training')
    ax.plot(np.log(E_test),label='test')
    ax.set_xlabel('iterations')
    ax.set_ylabel('log(E)')
    fig.suptitle(title)
    ax.legend()
    return fig

def save_figures(fig,path):
    for i,figure in enumerate(fig):
        figure.savefig(path+'{}'.format(i))







#These are the learning algorithms for the exercise 

def gradient_descent(w1,X,t,X_test,t_test,e=0.05,runs=10000): #e=learning_rate
    T1 = time.time()
    w = w1.copy()
    E_training = []
    E_test = []
    for i in range(runs):
        w += -e*dE(w,X,t)
        E_training.append(E(w,X,t))
        E_test.append(E(w,X_test,t_test))
        if i%100==0:
            T2 = time.time()
            print(i,str(round(T2-T1,1))+'s','E_training='+str(round(E(w,X,t),3)))
    T2 = time.time()
    dt = T2-T1
    return w, dt,E_training,E_test

def momentum(w0,X,t,X_test,t_test,e=0.05,a=0.03,runs=10000): #e=learning_rate, a=momentum_strength, 
    T1 = time.time()
    w = w0.copy()
    dw_old = 0
    E_training = []
    E_test = []
    for i in range(runs):
        dw_new = -e*dE(w,X,t) + a*dw_old
        w += dw_new + dw_old
        E_training.append(E(w,X,t))
        E_test.append(E(w,X_test,t_test))
        dw_old = dw_new.copy()
        if i%100==0:
            T2 = time.time()
            print(i,str(round(T2-T1,1))+'s','E_training='+str(round(E(w,X,t),3)))
    T2 = time.time()
    dt = T2-T1
    return w, dt,E_training,E_test

def weight_decay(w0,X,t,X_test,t_test,e=0.05,a=0.03,k=0.01,runs=10000): #e=learning_rate, a=momentum_strength, k=weight_decay_factor
    T1 = time.time()
    w = w0.copy()
    dw_old = 0
    E_training = []
    E_test = []
    for i in range(runs):
        dw_new = -e*dE(w,X,t) + a*dw_old
        w += dw_new + dw_old
        E_training.append(E(w,X,t))
        E_test.append(E(w,X_test,t_test))
        dw_old = dw_new.copy()
        if i%100==0:
            T2 = time.time()
            print(i,str(round(T2-T1,1))+'s','E_training='+str(round(E(w,X,t),3)))
    T2 = time.time()
    dt = T2-T1
    return w, dt,E_training,E_test

def newton_method(w0,X,t,X_test,t_test,e,a,k,runs): #e=learning_rate, a=momentum_strength, k=weight_decay_factor
    T1 = time.time()
    w = w0.copy()
    E_training = []
    E_test = []
    for i in range(runs):
        w += -np.matmul(inv(H_weight_decay(w,X,t,k)),dE_weight_decay(w,X,t,k)) 
        E_training.append(E(w,X,t))
        E_test.append(E(w,X_test,t_test))
        if i%1==0: 
            T2 = time.time()
            print(i,str(round(T2-T1,1))+'s','E_training='+str(round(E(w,X,t),3)))
    T2 = time.time()
    dt = T2-T1
    return w, dt,E_training,E_test

def line_search(w1,X,t,X_test,t_test,runs=300): 
    T1 = time.time()
    w = w1.copy()
    E_training = []
    E_test = []
    for i in range(runs):
        d = -dE_line_search(w)
        g = optimize.line_search(E_line_search,dE_line_search,w,d)[0]
        w += g*d
        E_training.append(E(w,X,t))
        E_test.append(E(w,X_test,t_test))
        if i%10==0:
            T2 = time.time()
            print(i,str(round(T2-T1,1))+'s','E_training='+str(round(E(w,X,t),3)))
    T2 = time.time()
    dt = T2-T1
    return w, dt,E_training,E_test

def conjugate_gradient_descent(w1,X,t,X_test,t_test,runs=200):
    T1 = time.time()
    w_old = w1.copy() #initialize both on w_old and w(=w_new) on the w1 that you start with 
    w = w1.copy()
    d = 0 #initialize d as 0 such that d is the gradient descent the first time it is used 
    E_training = []
    E_test = []
    for i in range(runs):
        b = dot(dE(w,X,t)-dE(w_old,X,t),dE(w,X,t))/np.dot(dE(w_old,X,t),dE(w_old,X,t))
        d = -dE(w,X,t) + b*d
        g = optimize.line_search(E_line_search,dE_line_search,w,d)[0]
        w_old = w.copy()
        w += g*d
        E_training.append(E(w,X,t))
        E_test.append(E(w,X_test,t_test))
        if i%10==0:
            T2 = time.time()
            print(i,str(round(T2-T1,1))+'s','E_training='+str(round(E(w,X,t),3)))
    T2 = time.time()
    dt = T2-T1
    return w, dt,E_training,E_test

def stochastic_gradient_descent(w1,X1,t1,X_test,t_test,e=0.01,div_factor=100,runs=5000): 
    T1 = time.time()
    w = w1.copy()
    Xi = np.array(np.array_split(X1,div_factor),dtype=object)
    ti = np.array(np.array_split(t1,div_factor),dtype=object)
    E_training = []
    E_test = []
    for i in range(runs):
        j = np.random.randint(0,div_factor)
        w += -e*dE(w,Xi[j],ti[j]) 
        E_training.append(E(w,X,t))
        E_test.append(E(w,X_test,t_test))
        if i%100==0:
            T2 = time.time()
            print(i,str(round(T2-T1,1))+'s','E_training='+str(round(E(w,X,t),3))) 
    T2 = time.time()
    dt = T2-T1
    return w, dt,E_training,E_test





#This is executing the script 
w0 = np.random.normal(0,1,size=X.shape[1]) #initializing a random model w0 
fig = []

w_gradient_descent,t_gradient_descent,E_training_gradient_descent, E_test_gradient_descent = gradient_descent(w0,X,t,X_test,t_test) #learning the model via gradient descent 
results_gradient_descent = [[E(w_gradient_descent,X,t),wrong_patterns(w_gradient_descent,X,t)],[E(w_gradient_descent,X_test,t_test),wrong_patterns(w_gradient_descent,X_test,t_test)]]
fig_gradient_descent = E_figures(E_training_gradient_descent,E_test_gradient_descent,'gradient descent')
fig.append(fig_gradient_descent)

w_momentum, t_momentum,E_training_momentum, E_test_momentum = momentum(w0,X,t,X_test,t_test) #learning the model via momentum 
results_momentum = [[E(w_momentum,X,t),wrong_patterns(w_momentum,X,t)],[E(w_momentum,X_test,t_test),wrong_patterns(w_momentum,X_test,t_test)]]
fig_momentum = E_figures(E_training_momentum, E_test_momentum,'momentum')
fig.append(fig_momentum)

w_weight_decay, t_weight_decay,E_training_weight_decay, E_test_weight_decay = weight_decay(w0,X,t,X_test,t_test) #learning the model via weight decay 
results_weight_decay = [[E(w_weight_decay,X,t),wrong_patterns(w_weight_decay,X,t)],[E(w_weight_decay,X_test,t_test),wrong_patterns(w_weight_decay,X_test,t_test)]]
fig_weight_decay = E_figures(E_training_weight_decay, E_test_weight_decay,'weight decay')
fig.append(fig_weight_decay)

#w_newton_method, t_newton_method,E_training_newton_method, E_test_newton_method = newton_method(w0,X,t,X_test,t_test,e=0.05,a=0.03,k=0.01,runs=10) #learning the model via weight decay 
#results_newton_method = [[E(w_newton_method,X,t),wrong_patterns(w_newton_method,X,t)],[E(w_newton_method,X_test,t_test),wrong_patterns(w_newton_method,X_test,t_test)]]
#fig_newton_method = E_figures(E_training_newton_method, E_test_newton_method,'newton method')
#fig.append(fig_newton_method)

w_line_search, t_line_search,E_training_line_search, E_test_line_search = line_search(w0,X,t,X_test,t_test) #learning the model via weight decay 
results_line_search = [[E(w_line_search,X,t),wrong_patterns(w_line_search,X,t)],[E(w_line_search,X_test,t_test),wrong_patterns(w_line_search,X_test,t_test)]]
fig_line_search = E_figures(E_training_line_search, E_test_line_search,'line search')
fig.append(fig_line_search)

#w_conjugate_gradient_descent, t_conjugate_gradient_descent,E_training_conjugate_gradient_descent, E_test_conjugate_gradient_descent = conjugate_gradient_descent(w0,X,t,X_test,t_test) #learning the model via weight decay 
#results_conjugate_gradient_descent = [[E(w_conjugate_gradient_descent,X,t),wrong_patterns(w_conjugate_gradient_descent,X,t)],[E(w_conjugate_gradient_descent,X_test,t_test),wrong_patterns(w_conjugate_gradient_descent,X_test,t_test)]]
#fig_conjugate_gradient_descent = E_figures(E_training_conjugate_gradient_descent, E_test_conjugate_gradient_descent,'conjugate gradient descent')
#fig.append(fig_conjugate_gradient_descent)

w_stochastic_gradient_descent, t_stochastic_gradient_descent,E_training_stochastic_gradient_descent, E_test_stochastic_gradient_descent = stochastic_gradient_descent(w0,X,t,X_test,t_test) #learning the model via weight decay 
results_stochastic_gradient_descent = [[E(w_stochastic_gradient_descent,X,t),wrong_patterns(w_stochastic_gradient_descent,X,t)],[E(w_stochastic_gradient_descent,X_test,t_test),wrong_patterns(w_stochastic_gradient_descent,X_test,t_test)]]
fig_stochastic_gradient_descent = E_figures(E_training_stochastic_gradient_descent, E_test_stochastic_gradient_descent,'stochastic gradient descent')
fig.append(fig_stochastic_gradient_descent)


#show results
print('gradient_descent')
show_results(results_gradient_descent,t_gradient_descent)
print('momentum')
show_results(results_momentum,t_momentum)
print('weight decay')
show_results(results_weight_decay,t_weight_decay)
print('line search')
show_results(results_line_search,t_line_search)
print('stochastic gradient descent')
show_results(results_stochastic_gradient_descent,t_stochastic_gradient_descent)

#save_figures(fig,'/Users/samsuidman/Desktop/neurophysics/machine_learning/machine_learning_wc4_E_decay_') #save the figures

#visualize(w_grad_descent,X,1328) #visualizing a learned model for a test handwritten digit 


