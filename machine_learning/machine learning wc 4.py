import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.io import loadmat


file = loadmat('/Users/samsuidman/Desktop/neurophysics/machine_learning/mnistAll.mat') 

train_images = file['mnist']['train_images'][0][0]
test_images = file['mnist']['test_images'][0][0]
train_labels = file['mnist']['train_labels'][0][0].transpose()[0]
test_labels = file['mnist']['test_labels'][0][0].transpose()[0]


indices_3 = np.where(train_labels==3)[0]
indices_7 = np.where(train_labels==7)[0]
X3 = train_images[:,:,indices_3]
X7 = train_images[:,:,indices_7]
n3 = np.size(X3,2)
n7 = np.size(X7,2)
X3 = np.reshape(X3,[784,n3])
X7 = np.reshape(X7,[784,n7])
X3 = X3/np.max((np.max(np.concatenate(X3)),np.max(np.concatenate(X7))))
X7 = X7/np.max((np.max(np.concatenate(X3)),np.max(np.concatenate(X7))))
X3 = np.insert(X3,0,1,axis=1)
X7 = np.insert(X3,0,1,axis=1)
X = [X3,X7]
t = [np.zeros([1,n3])[0],np.ones([1,n7])[0]]





