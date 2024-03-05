import numpy as np
import matplotlib.pyplot as plt

#4b

def err(N,P_training): #making a function so N and P_training can be adapted
    #set some constant variables
    P_test = 10000 #amount of test sets 
    n_runs = 100 #amount of runs a model w to be learned
    e = 0.01 #learning rate
    n_w = 20 #amount of models that need to be learned 

    #create training set
    x_training = np.random.choice([-1,1],size=(P_training,N)) #create P training sets

    #create teacher values 
    w_teacher = np.random.normal(0,1,size=N) #create a (random) model that we need to find later 
    y_teacher = np.sign(np.matmul(x_training,w_teacher)) #with this model certain values of y can be created

    #create models w that need be learned
    w = np.random.normal(0,1,size=[n_w,N]) #Now random vectors of w are created. These vectors are different and they all need to learn the model. This way you can see different w's that will be learned, but they still converge. Eventually with these converged w's the x_training will all get to y_training so gp(f)=1. With this we can look at g(f) and the generalization bound. 

    #training the models w on the training set 
    convergence = np.array([])
    for j in range(n_w): #to training all the different w's
        p = 0
        for i in range(n_runs):
            y_i = np.heaviside(-np.dot(x_training,w[j])*y_teacher,0) #dot products are calculated for each pattern with a w (the one that is trained at the moment). The ones that are negative need to be updated, so they get a 1. This y_i has 1000 elements, a single 0 or 1 for each pattern. The 1 is if a pattern is not yet predicted correct.  
            dw = (x_training.transpose()*y_teacher*y_i).transpose().sum(axis=0) #add all patterns times y values (so you get the x_mu in the slides) (for each element in the vector) together, because they all need to be added to the w (with a learning factor in front). 
            w[j] += e*dw #the w that is trained at the moment needs to be updated by all the patterns that have a negative dot product. The learning factor determines how fast this happens. 
            if np.sum(y_i)==0 and p==0:
                convergence = np.append(convergence,i+1)
                p = 1

    #checking if the models w indeed give the right result for the teacher set
    check_array = np.array([]) #here you can check if all the w that you have trained really get the y value that they have been trained to. 
    for j in range(n_w): #look at all w's
        check_value = np.all(np.sign(np.dot(x_training,w[j])) == y_teacher) #check of all y are predicted correctly by a model w 
        check_array = np.append(check_array,check_value) #make an array that needs to return only 1's (so True values)

    #creating test set and test values 
    x_test = np.random.choice([-1,1],size=(P_test,N)) 
    y_teacher_test = np.sign(np.matmul(x_test,w_teacher)) 

    #calculate the errors of the sets 
    errors = np.array([])
    for j in range(n_w):
        y = np.sign(np.dot(x_test,w[j])) #calculating the y value on the test sets for a certain learned w 
        error = np.average(y!=y_teacher_test) #calculating how much of the predicted y's are wrong compared to the created set. The average of the amount of wrong sets is the error
        errors = np.append(errors,error) #add the error to the error list for all models w
    
    #average the errors 
    estimated_error = np.average(errors) #calculate the average error over all the models w

    return estimated_error #return the error 



def calculate(N,P,delta):
    eta = np.sqrt(8/P*np.log(4/delta*((2*np.e*P)/N)**N))
    return eta


errors_numerical_estimates = np.array([]) #the trained errors for different P_training
error_generalization_bound = np.array([]) #the calculated errors for different P_training
for p in [10,50,100,500,1000]: 
    value1 = err(10,p) #training function
    value2 = calculate(10,p,0.01) #calculating function
    errors_numerical_estimates = np.append(errors_numerical_estimates,value1)
    error_generalization_bound = np.append(error_generalization_bound,value2)



