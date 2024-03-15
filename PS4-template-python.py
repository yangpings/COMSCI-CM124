#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np 
import pandas as pd

import matplotlib.pyplot as plt
from scipy.special import logsumexp
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score


# # read data

# In[2]:


X_df = pd.read_csv("./data/Q3/mixture1.geno", header = None)
Xnew_df = pd.read_csv("./data/Q3/mixture2.geno", header = None)
F_df = pd.read_csv("./data/Q3/mixture1.freq", header = None)
Z_df = pd.read_csv("./data/Q3/mixture1.ganc", header = None)


# # implementation M

# In[3]:


#X: N by M data matrix
#gamma: P(Z_i=k | X_i=x_i, theta) N by K
def M_step(X, gamma): 
    N, M = X.shape
    K = gamma.shape[1]
    
    ######### TODO 3a: modify the following to have meaning updates #########
    pis = None #Proportion vector pi: length K
    F = matrix(0, M, K) #Frequency matrix F: M by K
    
    
    
    
    
    ######### end of modification #########     
    return ({"pis": pis, "F": F})


# # implementation E

# In[4]:


#X: N by M data matrix
#params: a dictionary with two parameters returned from M_step
def E_step(X, params, thr = 10**(-8)):
    F = params["F"] #Frequency matrix F: M by K
    pis = params["pis"] #Proportion vector pi: length K

    N, M = X.shape
    K = F.shape[1] 

    ######### TODO 3b: modify the following to have meaning updates #########
    #calculate weighted_log_prob: log(P(X_i=x_i | Z_i=k, theta) * P(Z_i=k | theta))
    #calcualte log_prob_sample: log P(Xi=x_i | theta) length N vector. Hint: use logsumexp function
    #calcualte log_prob_data: log P(X_1:n=x_1:n | theta) scalar
    #calculate log_gammas: log P(Z_i=k | X_i=x_i, theta) N by K
    weighted_log_prob = np.zeros((N,K))
    
    
    
    
    

    ######### end of modification #########
    return log_gammas, log_prob_data


# # implementation EM 

# In[5]:


def EM(X, K = 2, max_iter = 100, tol = 10**(-4), n_init = 3, debug = False):
    
    N, M = X.shape 
    res = {}
    best_log_prob_data = -np.inf
    converged  = False 

    #loop through different random starting points
    for init in range(1, 1+n_init, 1):
        np.random.seed(init)
        if(debug):
            print(f"starting EM on random initialization: {init} out of {n_init}")
        
        ######### TODO 3c: modify the following to have the full EM updates #########
        
        #initialize soft assignment 
        gammas = None
        
    
        
        log_prob_data  = -np.inf
        for n_iter in range(1, 1+max_iter, 1):
            prev_log_prob_data = log_prob_data
        
        
        
        
        
        
        
            ######### convergence check #########
            change = (log_prob_data - prev_log_prob_data)/N
            if abs(change) < tol:
                if(debug):
                    print(f"random initialization {init} converged at iteration {n_iter}")
                    print("")
                converged = True
                break


        
        
        
        ######### update on the best initialization #########
        best_init = NULL
        
        
        
        
        
    ######### end of modification #########
    res["converged"] = converged     
    res["best_init"] = best_init
    return(res)

