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

def E_line_search(w): #Return the error for w,X and the actual value 3,7 (so actual 1,0)
    y0 = 1/(1+np.exp(-np.dot(X0,w)))
    Ew = -sum(t0*log(y0)+(1-t0)*log(1-y0))/len(X0)
    return Ew

def dE_line_search(w): #gradient for w,X,t
    y0 = 1/(1+np.exp(-np.dot(X0,w)))
    dE = matmul((y0-t0),X0)/len(X0)
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





#These are the learning algorithms for the exercise 

def grad_descent(w0,X,t,e,runs): #e=learning_rate
    T1 = time.time()
    w = w0.copy()
    for i in range(runs):
        w += -e*dE(w,X,t)
        if i%100==0:
            T2 = time.time()
            print(i,str(round(T2-T1,1))+'s','E='+str(round(E(w,X,t),3)))
    return w

def momentum(w0,X,t,e,a,runs): #e=learning_rate, a=momentum_strength, 
    T1 = time.time()
    w = w0.copy()
    dw_old = 0
    for i in range(runs):
        dw_new = -e*dE(w,X,t) + a*dw_old
        w += dw_new + dw_old
        dw_old = dw_new.copy()
        if i%100==0:
            T2 = time.time()
            print(i,str(round(T2-T1,1))+'s','E='+str(round(E(w,X,t),3)))
    return w

def weight_decay(w0,X,t,e,a,k,runs): #e=learning_rate, a=momentum_strength, k=weight_decay_factor
    T1 = time.time()
    w = w0.copy()
    dw_old = 0
    for i in range(runs):
        dw_new = -e*dE(w,X,t) + a*dw_old
        w += dw_new + dw_old
        dw_old = dw_new.copy()
        if i%100==0:
            T2 = time.time()
            print(i,str(round(T2-T1,1))+'s','E='+str(round(E(w,X,t),3)))
    return w

def newton_method(w0,X,t,e,a,k,runs): #e=learning_rate, a=momentum_strength, k=weight_decay_factor
    T1 = time.time()
    w = w0.copy()
    for i in range(runs):
        w += -np.matmul(inv(H_weight_decay(w,X,t,k)),dE_weight_decay(w,X,t,k))
        if i%1==0: 
            T2 = time.time()
            print(i,str(round(T2-T1,1))+'s','E='+str(round(E(w,X,t),3)))
    return w

def line_search(w,X,t,runs):
    T1 = time.time()
    w = w0.copy()
    for i in range(runs):
        d = -dE_line_search(w)
        g = optimize.line_search(E_line_search,dE_line_search,w0,d)[0]
        w += g*d
        if i%10==0:
            T2 = time.time()
            print(i,str(round(T2-T1,1))+'s','E='+str(round(E(w,X,t),3)))
    return w

def conjugate_gradient_descent(w,X,t,runs):
    T1 = time.time()
    w = w0.copy()
    w_old = w.copy()
    d = -dE(w,X,t)
    for i in range(runs):
        b = dot(dE(w,X,t)-dE(w_old,X,t),dE(w,X,t))/np.dot(dE(w_old,X,t),dE(w_old,X,t))
        d = -dE(w,X,t) + b*d
        errors = np.array([E(w+i*d,X,t) for i in np.linspace(0.05,0.9,30)]) #create array with all errors for a certain gammas. 
        min_error = np.where(errors == np.min(errors))[0][0] #find the index of the minimal error. 
        g = np.linspace(0.05,0.9,30)[min_error] 
        w_old = w.copy()
        w += g*d
        if i%10==0:
            T2 = time.time()
            print(i,str(round(T2-T1,1))+'s','E='+str(round(E(w,X,t),3)))
            print(g)
    return w

def stochastic_gradient_descent(w0,X0,t0,e,div_factor,runs): 
    T1 = time.time()
    w = w0.copy()
    Xi = np.array(np.array_split(X,div_factor),dtype=object)
    ti = np.array(np.array_split(t,div_factor),dtype=object)
    for i in range(runs):
        j = np.random.randint(0,div_factor)
        w += -e*dE(w,Xi[j],ti[j]) #waarschijnlijk gaat het in deze stap fout, omdat de som van al die gradient descent termen heel groot is. 
        if i%100==0:
            T2 = time.time()
            print(i,str(round(T2-T1,1))+'s','E='+str(round(E(w,X,t),3))) #Ook kijken naar de errorfunctie 
    return w





#This is executing the script 
w0 = np.random.normal(0,1,size=X.shape[1]) #initializing a random model w0 

w_grad_descent = grad_descent(w0,X,t,e=0.05,runs=10000) #learning the model via gradient descent 
results_grad_descent = [[E(w_grad_descent,X,t),wrong_patterns(w_grad_descent,X,t)],[E(w_grad_descent,X_test,t_test),wrong_patterns(w_grad_descent,X_test,t_test)]]

w_momentum = momentum(w0,X,t,e=0.05,a=0.03,runs=10000) #learning the model via momentum 
results_momentum = [[E(w_momentum,X,t),wrong_patterns(w_momentum,X,t)],[E(w_momentum,X_test,t_test),wrong_patterns(w_momentum,X_test,t_test)]]

w_weight_decay = weight_decay(w0,X,t,e=0.05,a=0.03,k=0.01,runs=10000) #learning the model via weight decay 
results_weight_decay = [[E(w_weight_decay,X,t),wrong_patterns(w_weight_decay,X,t)],[E(w_weight_decay,X_test,t_test),wrong_patterns(w_weight_decay,X_test,t_test)]]

w_newton_method = newton_method(w0,X,t,e=0.05,a=0.03,k=0.01,runs=10) #learning the model via weight decay 
results_newton_method = [[E(w_newton_method,X,t),wrong_patterns(w_newton_method,X,t)],[E(w_newton_method,X_test,t_test),wrong_patterns(w_newton_method,X_test,t_test)]]

X0 = X.copy()
t0 = t.copy()
w_line_search = line_search(w0,X,t,runs=200) #learning the model via weight decay 
results_line_search = [[E(w_line_search,X,t),wrong_patterns(w_line_search,X,t)],[E(w_line_search,X_test,t_test),wrong_patterns(w_line_search,X_test,t_test)]]

w_conjugate_gradient_descent = conjugate_gradient_descent(w0,X,t,runs=200) #learning the model via weight decay 
results_conjuage_gradient_descent = [[E(w_conjugate_gradient_descent,X,t),wrong_patterns(w_conjugate_gradient_descent,X,t)],[E(w_conjugate_gradient_descent,X_test,t_test),wrong_patterns(w_conjugate_gradient_descent,X_test,t_test)]]

w_stochastic_gradient_descent = stochastic_gradient_descent(w0,X,t,e=0.01,div_factor=100,runs=5000) #learning the model via weight decay 
results_stochastic_gradient_descent = [[E(w_stochastic_gradient_descent,X,t),wrong_patterns(w_stochastic_gradient_descent,X,t)],[E(w_stochastic_gradient_descent,X_test,t_test),wrong_patterns(w_stochastic_gradient_descent,X_test,t_test)]]


#visualizing somehting 
visualize(w_grad_descent,X,1328) #visualizing a learned model for a test handwritten digit 


