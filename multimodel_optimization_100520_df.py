#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


#%%


ensemble_size=10                                         # size of ensemble at any time
num_replace_different=2                                  # number of metropolis models to add because they are different from previously sampled
num_replace_good=1                                       # number of metropolis models to add because they are close to the retained models in the ensemble
metropolis_rounds=20                                     # number of rounds to try metropolis replacement
Pnum=ensemble_size*(metropolis_rounds+1)                 # number of parameter sets to consider (size of sampled world)
type(Pnum)


#%%


Params = ('v', 'D', 'Lambda', 'R','xsampled')   
type(Params)                
 # xsampled is a flag to indicate which samples have been considered already
Plims=pd.DataFrame([[1,5],[4,9],[5,23],[23,27],[0,0]])           # parameter limits
Plims.index = Params
Plims.columns = ('minval','maxval')
Plims

#%%

Params = ('v', 'D', 'Lambda', 'R','xsampled','xID')
Pspace = pd.DataFrame(np.random.rand(Pnum,Plims.shape[0]+1),columns=Params)         # fill the sampled world with random parameter sets - may want to use Latin Hypercube Sampling here
Pspace = Pspace * (Plims.maxval-Plims.minval) + Plims.minval                       # convert random numbers to parameter values within ranges

#Explore previous line for generating breakthrough curve
#check by plotting the Pspace for distribution

Pspace.xsampled=0                                                                  # set sampled flag to zero initially
Pspace.xID=np.arange(Pnum) 


#%%
#New stuff



jump=10                                             #pulls ten sets of variables
intervals=Pspace[Pspace.xsampled==0].iloc[0:jump]   #Find where xsampled of Pspace equals 0, and take a "jump" amount of sets
Pspace.xsampled[0:jump]=1                           #set the "jump" amount of sets to have xsampled=1, to not sample again


#%%
#Next steps:

#run the model in a way to feed the entirety of intervals and get out a breakthrough curve for every row
#try using "intervals" in C equation to feed df all together


# In[12]:


Pspace.head(20) #check that the first set of "jump" have now been flagged as 1 and will not be resampled on the next run


# In[ ]:




