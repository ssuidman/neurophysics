import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import os.path as path

#Loading in the MNIST data
file = loadmat('/Users/samsuidman/Desktop/neurophysics/machine_learning/mnistAll.mat') 

#Preparing the data
train_images = file['mnist']['train_images'][0][0]
X = np.transpose(train_images.reshape(train_images.shape[0]*train_images.shape[1],train_images.shape[2]))
X = X/np.max(X)

N = X.shape[0]
d = X.shape[1]
K = 10

def get_k(X,mu): #hier komt een lijst van allemaal dezelfde k uit, waarschijnlijk door bepaalde initialisatie. Later kijken of dit fout is. 
    k_matrix = np.matmul(np.log(np.where(mu==0,1,mu)),X.transpose()) + np.matmul(np.log(1-np.where(mu==1,0,mu)),(1-X).transpose()) #Here mu can be 0 since sum of all first digits of image is zero because it is background. In this case x*log(mu) = 0*log(0) = 0. That is why mu should be set to 1 inside the log, so log(1)=0. 
    k = np.argmax(k_matrix,axis=0)
    return k 

def get_pi(k,K,N):
    pi = np.array([len(np.where(k==i)[0]) for i in range(K)])/N
    return pi

def get_mu(X,pi,k,d,K,N):
    mu = np.array([np.sum(X[np.where(k==i)],axis=0)/(pi[i]*N) if pi[i]!=0 else np.zeros(d) for i in range(K)])
    return mu 

def clustering(X,d,K,N,runs):
    k = np.random.randint(10,size=d)
    for i in range(runs):
        print(i)
        pi = get_pi(k,K,N)
        mu = get_mu(X,pi,k,d,K,N)
        k = get_k(X,mu)
    return k, pi, mu

def figures(K,n=10): #gives for each cluster a figure with nxn pictures
    fig = []
    for j in range(K):
        print(j)
        index = np.where(k==j)
        cluster = X[index]
        fig0, ax = plt.subplots(ncols=n,nrows=n)
        for i in range(n):
            for j in range(n):
                ax[i][j].imshow(cluster[10*i+j].reshape(28,28))
                ax[i][j].xaxis.set_visible(False)
                ax[i][j].yaxis.set_visible(False)
        fig.append(fig0)
    return fig

def save_figures(fig,path):
    for i,figure in enumerate(fig):
        figure.savefig(path+'{}'.format(i))
        
k,pi,mu = clustering(X,d,K,N,100) #get the variables k,pi,mu
fig = figures(K,n=10) #get some figures to visualize the clusters
save_figures(fig,'/Users/samsuidman/Desktop/neurophysics/machine_learning/machine_learning_wc7_cluster_') #save the clusters

#Total running time is about 3 minutes 
