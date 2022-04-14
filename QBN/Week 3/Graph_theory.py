import numpy as np
import matplotlib.pyplot as plt

# problem 2a
n = 30
m = 100
A = np.zeros(n**2)
A[:m] = 1
np.random.shuffle(A)
A = np.reshape(A,[n,n])
A_sum = np.sum(A,axis=0)
C = A_sum/n
plt.hist(C,bins=10)
print(A)

#problem 2b
n = 30
A = np.random.choice([0,1],size=[n,n],p=[770/870,100/870])
A_sum = np.sum(A,axis=0)
C = A_sum/n
plt.hist(C,bins=10)
# print(A)

#problem 2c
n = 14
k = 4
A = np.zeros([n,n])
for i in range(n):
    for j in range(n):
        if i!=j and (np.abs(i-j)<=k or np.abs(i-j)>=np.abs(n-k)):
            A[i][j] = 1
print(A)