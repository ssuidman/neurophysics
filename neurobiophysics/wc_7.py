def get_musclefit(target,P,exponent=1.1,lamda=0.5,stdev=0):
    alpha = 0.02     # Gradient descent learning rate 
    # Gradient descent algorithm on the cost finds the optimal u_i within maximum 5000 iterations
    #
    #  first, parameter initialization:
    n = 1                                             # iteration counter
    u = np.zeros(P.shape[1])                          # initialize the optimal ui to  zeros
    y = np.matmul(P,u)
    error = target - y 
    cost = np.dot(error,error) + lamda * np.sum(u**exponent) # the cost function (no noise)
    du = np.array([1,1])
    costgrad = np.zeros(P.shape[1]) # initialize the cost gradient to zero
    
    while (np.sqrt(np.dot(du,du))>0.000001 and n<5000): # stop criteria
        # costgrad = -2*P' * error +lamda * exponent * u.^(exponent-1);     
        # this is the cost gradient with respect to ui without signal-dependent noise
        n += 1
        for i in range(len(u)):
            costgrad[i] = -2*np.dot(error,P[:,i]) + 2*u[i]*stdev**2*np.dot(P[:,i],P[:,i]) + lamda*exponent*u[i]**(exponent-1)    
            # gradient is derived from equation 10.9, with signal-dependent noise included
        unew = u - alpha * costgrad
        unew = np.where(unew>0, unew, 0)        # enforce non-negativity constraint
        du = unew - u
        u = unew.copy()                                  # new estimate for the activations
        y = np.matmul(P,u)
        error = target - y                                  # update the error.
        cost = np.dot(error,error) + lamda * np.sum(u**exponent)
    return u,error

Pa = np.array([9, 49, 100, 145, 190, 230, 270, 340])/180*np.pi # P in radians
P=np.array([np.cos(Pa), np.sin(Pa)]) # 2D unit vectors
target = np.array([1,0])
u, error = get_musclefit(target,P)


for i in range(P.shape[1]):
    plt.plot([0,u[i]*P[0,i]],[0,u[i]*P[1,i]])
plt.plot([0,np.matmul(P,u)[0]],[0,np.matmul(P,u)[1]],label='sum force')
plt.plot([0,target[0]],[0,target[1]],label='target')
plt.legend()
