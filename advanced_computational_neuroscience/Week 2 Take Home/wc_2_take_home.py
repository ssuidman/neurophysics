import numpy as np
import matplotlib.pyplot as plt
from skimage import io, color

# DATA = 'pixels'
# DATA = 'gratings'
DATA = 'natural'

# METHOD = 'pca'
# METHOD = 'sparse'


# Generate data
if DATA == 'pixels':
    n1 = 16 # image size is n1 * n1
    N = 1500 # number of images
    n = n1**2
    x = np.zeros([N,n1,n1])
    y = np.random.rand(N,n)
    x = -np.log(y)
    x = np.reshape(x,[N,n1,n1])
    x /= np.max(x)

if DATA == 'gratings':
    n1 = 16 # image size is n1 * n1
    N = 1000 # number of images
    n = n1**2
    x = np.zeros([N,n1,n1])
    n2 = 4 # number of grating components
    k = 2*np.pi*np.array([1/2,1/4,1/8,1/16]) # grating wave number
    for i in range(int(N/2)):
        y = np.random.rand(n2)
        a = -np.log(y)
        t = np.outer(k,np.arange(1,n1+1))
        u = np.matmul(a,np.sin(t))
        x[i] = np.outer(np.ones(n1),u)
        x[int(N/2)+i] = np.outer(u,np.ones(n1))
    x /= np.max(x)

if DATA == 'natural':
    img = io.imread('/Users/samsuidman/Desktop/neurophysics/advanced_computational_neuroscience/Week 2 Take Home/sampleMerry_0011_Lasalle.jpeg')
    g = color.rgb2gray(img)
    plt.imshow(g,cmap='gray')
#     plt.savefig('test.png')
    a1,a2 = g.shape
    k = 0
    n1 = 16 
    n = n1**2
    N = int(a1*a2/n)
    x = np.zeros([N,n1,n1])
    for i in range(int(a1/n1)): 
        for j in range(int(a2/n1)): 
            x[k] = g[i*n1:(i+1)*n1,j*n1:(j+1)*n1]
            k += 1
    x /=265

I = x.reshape(N,n) # reshape X
I -= np.mean(I,axis=0) # substract mean, such that sum of images for each pixel is zero

C0 = np.matmul(np.transpose(I),I)/N
w,v = np.linalg.eig(C0) # w[i] eigenvalues with orthonormal eigenvectors v[:,i]
K = 25
i_PCA = np.argsort(w)[-K:][::-1]
phi = v[:,i_PCA] # with phi[:,i] the i'th principle component

a = np.matmul(I,np.transpose(v)) # to get the n weights for all N images 
I_check = np.matmul(a,v) # this way you get back the images via eigenvectors and weights 
np.sum(I_check-I) # this is a test to check that you indeed get back the images with an error around 10^-17 

fig,ax = plt.subplots(nrows=5,ncols=5)
for i in range(5):
    for j in range(5):
        ax[i][j].imshow(v[:,i*3+j].reshape(n1,n1),cmap='gray')
fig.colorbar(plt.cm.ScalarMappable(cmap='gray'),ax=ax)









#####################################################################
################ fout 1 ##############################################
#####################################################################
eta = 0.0001 # learning rate
lamda = 0.14 # penalty

En = np.sum(1/2*(I-np.matmul(a,phi))**2,axis=1) + lamda*np.sum(np.abs(a),axis=1)# calculate error En[pictures]
E = np.mean(En) # mean error between pictures 

b = np.matmul(I,np.transpose(phi)) # b[pictures,features] = sum( I[pictures,pixels]*phi[features,pixels]^T ) 
C = np.matmul(phi,np.transpose(phi)) # C[features,features] 
dE_a = -b + np.matmul(a,C) + lamda*np.sign(a) # dE_a[pictures,features] 
dE_phi = -np.matmul(np.transpose(a),I-np.matmul(a,phi)) # dE_phi[features,pixels]    =   sum( [pictures,features]^T*[pictures,pixels] ) 

for i in range(500): 
    dE_phi = -np.matmul(np.transpose(a),I-np.matmul(a,phi)) # dE_phi[features,pixels] 
#     phi += -eta*dE_phi 
    
    b = np.matmul(I,np.transpose(phi)) # b[pictures,features] = sum( I[pictures,pixels]*phi[features,pixels]^T ) 
    C = np.matmul(phi,np.transpose(phi)) # C[features,features] 
    a = np.where(np.abs(b-np.matmul(a,C)+a*np.diag(C))<=lamda,0,b-lamda*np.sign(b)) # a[pictures,features]
    En = np.sum(1/2*(I-np.matmul(a,phi))**2,axis=1) + lamda*np.sum(np.abs(a),axis=1)# calculate error En[pictures] 
    E = np.mean(En) # total error
    if i in 100*np.arange(10):
        print(i,E)
print('done')










###########################################################################################
####################################### fout 2 ############################################
###########################################################################################
I = x.reshape(N,n) # reshape X
I -= np.mean(I,axis=0) # substract mean, such that sum of images for each pixel is zero

C0 = np.matmul(np.transpose(I),I)/N
w,v = np.linalg.eig(C0) # w[i] eigenvalues with orthonormal eigenvectors v[:,i]
K = 25
i_PCA = np.argsort(w)[-K:][::-1]
phi = v[:,i_PCA] # with phi[:,i] the i'th principle component
# a = np.random.rand(N,n) # a[pictures,features]
a = np.random.normal(0,0.17,size=[N,n]) # a[pictures,features]

eta = 0.0001 # learning rate
lamda = 0.14 # penalty

En = np.sum(1/2*(I-np.matmul(a,phi))**2,axis=1) + lamda*np.sum(np.abs(a),axis=1)# calculate error En[pictures]
E = np.mean(En) # mean error between pictures 

# for i in range(500):     
#     b = np.matmul(I,np.transpose(phi)) # b[pictures,features] = sum( I[pictures,pixels]*phi[features,pixels]^T ) 
#     C = np.matmul(phi,np.transpose(phi)) # C[features,features] 
#     a = np.where(np.abs(b-np.matmul(a,C)+a*np.diag(C))<=lamda,0,b-lamda*np.sign(b)) # a[pictures,features]
#     En = np.sum(1/2*(I-np.matmul(a,phi))**2,axis=1) + lamda*np.sum(np.abs(a),axis=1)# calculate error En[pictures] 
#     E = np.mean(En) # total error
#     if i in 100*np.arange(10):
#         print(i,E)
# print('done')








###########################################################################################
####################################### fout 3 ############################################
############################## overleg dit met wc assistenten #############################
##### Het is gek dat je <=lamda doet, omdat je dan het beste gewoon lamda=0 kan kiezen ####
#### Je krijgt dan bijna alleen maar nullen voor, omdat je dat door de sparsity kiest. ####
#### Voor minimale E is a 100% 0 als je verschillende lambda kiest, dit is wel gek. #######
###########################################################################################

I = x.reshape(N,n) # reshape x to get I
I -= np.mean(I,axis=0) # I[pictures,pixels]. Substract mean, such that sum of images for each pixel is zero

C0 = np.matmul(np.transpose(I),I)/N # C0[pixels,pixels] 
w,v = np.linalg.eig(C0) # w[features] eigenvalues with orthonormal eigenvectors v[pixels,features] 
K = 25 
i_PCA = np.argsort(w)[-K:][::-1] 
phi = np.transpose(v[:,i_PCA]) # phi[features,pixels] 
b = np.matmul(I,np.transpose(phi)) # b[pictures,features] = sum( I[pictures,pixels]*phi[features,pixels]^T ) 
C = np.matmul(phi,np.transpose(phi)) # C[features,features] 
lamdas = np.array([0,0.000001,0.00001,0.0001,0.001,0.01,0.1])
E_array = np.zeros(len(lamdas))

for i,lamda in enumerate(lamdas):
# lamda = 0.014 # penalty
    a = np.where(np.abs(b)<=lamda,0,b-lamda*np.sign(b)) # a[pictures,features] 
    En = np.sum(1/2*(I-np.matmul(a,phi))**2,axis=1) + lamda*np.sum(np.abs(a),axis=1)# calculate error En[pictures]
    E = np.mean(En) # mean error between pictures 
    print(E)
    print( str(round((1-np.sum(np.where(a==0,0,1))/(N*K))*100,2)) + "%" )
    E_array[i] = E
print('Done!')
plt.plot(E_array)
