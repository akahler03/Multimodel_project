#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math
import matplotlib.pyplot as plt
import warnings
import pandas as pd
warnings.filterwarnings("ignore")
initial_c=1
x_obs=200
pulse_time=200
t_start=0               
t_end = 800
dt = 10
t_obs=np.arange(t_start,t_end,dt)  


# In[2]:


trials=5   #np.linspace(0,100,100)
vmin=1
vmax=3
Dmin=10
Dmax=100
Rmin=1
Rmax=2
decaymin=0
decaymax=.005


# In[3]:


v=np.random.uniform((vmin),(vmax),size=trials).round(4)
D=np.random.uniform((Dmin),(Dmax),size=trials).round(4)
R=np.random.uniform((Rmin),(Rmax),size=trials).round(4)
decay=np.random.uniform((decaymin),(decaymax),size=trials).round(4)
vtrue=1.6
Dtrue=55
Rtrue=1.3
decaytrue=0


# In[4]:


originals=pd.DataFrame(np.hstack((v[:,None], D[:,None],R[:,None],decay[:,None])))
originals.columns=['velocity','Dispersion','Retardation','decay']
originals.loc[ : , 'sampled'] = False
originals.loc[:,'mismatch']=np.nan


# In[5]:


originals

C=np.zeros((np.shape(t_obs)[0],np.shape(R)[0]))
countR=-1
for Rval in R:
    countR+=1
    C[:,countR]=((initial_c*pulse_time*x_obs*Rval**.5)/(2*(math.pi*D[countR]*(t_obs**3))**.5))*np.exp(-((Rval*x_obs-v[countR]*t_obs)**2)/(4*Rval*D[countR]*t_obs))*np.exp(-decay[countR]*t_obs)                          
C[np.where(C<0)]=0
C[np.where(np.isinf(C)==1)]=0
C[np.where(np.isnan(C)==1)]=0

C_true=np.zeros((np.shape(t_obs)[0]))
C_true[:]=((initial_c*pulse_time*x_obs*Rtrue**.5)/(2*(math.pi*Dtrue*(t_obs**3))**.5))*np.exp(-((Rtrue*x_obs-vtrue*t_obs)**2)/(4*Rtrue*Dtrue*t_obs))*np.exp(-decaytrue*t_obs)                          
C_true[np.where(C_true<0)]=0
C_true[np.where(np.isinf(C_true)==1)]=0
C_true[np.where(np.isnan(C_true)==1)]=0

err_val=np.random.normal(loc=0,scale=.03,size=(C_true.size))
C_true=C_true+err_val
# In[6]:


C_true=np.zeros((np.shape(t_obs)[0]))
C_true[:]=((initial_c*pulse_time*x_obs*Rtrue**.5)/(2*(math.pi*Dtrue*(t_obs**3))**.5))*np.exp(-((Rtrue*x_obs-vtrue*t_obs)**2)/(4*Rtrue*Dtrue*t_obs))*np.exp(-decaytrue*t_obs)                          
C_true[np.where(C_true<0)]=0
C_true[np.where(np.isinf(C_true)==1)]=0
C_true[np.where(np.isnan(C_true)==1)]=0

err_val=np.random.normal(loc=0,scale=.03,size=(C_true.size))
C_true=C_true+err_val


# In[7]:


originals.sampled


# In[ ]:




#FUNCTIONAL NEW CODE---SAVE

C=np.zeros((np.shape(t_obs)[0],trials))
count1=0
count2=10
nummetropolis=1    
for i in np.arange(nummetropolis):
    for j in np.arange(count1,count2):
        if originals.sampled[j]==False:
            v=originals.velocity[j]
            D=originals.Dispersion[j]
            R=originals.Retardation[j]
            decay=originals.decay[j]
            originals.sampled[j]==True
            C[:,j]=((initial_c*pulse_time*x_obs*R**.5)/(2*(math.pi*D*(t_obs**3))**.5))\
            *np.exp(-((R*x_obs-v*t_obs)**2)/(4*R*D*t_obs))*np.exp(-decay*t_obs) 
            
            #Insert something in this loop to change originals.sampled column to 1 for the range
            
            #Another "if" statement here to calculate mismatch, weigh against the next in range, and replace if better
            
            labelstr='run ' + str(j+1)
            plt.plot(t_obs,C[:,:],label=labelstr)
            plt.plot(t_obs,C_true,'.',label='C true w/error')
            plt.xlabel('Time (days)')
            plt.ylabel('Concentration') 
            plt.title('X_obs =  200m') 
            plt.show()       #Will plot one iteration at a time, cumulatively
# In[8]:


#CREATE NEW NEW CODE
C=np.zeros((np.shape(t_obs)[0],trials))


# In[ ]:





# In[9]:


#EDIT THE NEW CODE IN THIS CELL

C=np.zeros((np.shape(t_obs)[0],trials))
count1=0
count2=10
nummetropolis=1    
for i in np.arange(nummetropolis):
    for j in np.arange(count1,count2):
        if originals.sampled[j]==False:
            v=originals.velocity[j]
            D=originals.Dispersion[j]
            R=originals.Retardation[j]
            decay=originals.decay[j]
            originals.sampled[j]==True
            C[:,j]=((initial_c*pulse_time*x_obs*R**.5)/(2*(math.pi*D*(t_obs**3))**.5))            *np.exp(-((R*x_obs-v*t_obs)**2)/(4*R*D*t_obs))*np.exp(-decay*t_obs) 
            
            #Insert something in this loop to change originals.sampled column to 1 for the range
            
            #Another "if" statement here to calculate mismatch, weigh against the next in range, and replace if better
            
            labelstr='run ' + str(j+1)
            plt.plot(t_obs,C[:,:],label=labelstr)
            plt.plot(t_obs,C_true,'.',label='C true w/error')
            plt.xlabel('Time (days)')
            plt.ylabel('Concentration') 
            plt.title('X_obs =  200m') 
            #plt.show()       #Will plot one iteration at a time, cumulatively


# In[10]:


originals.sampled


# In[11]:


test=((initial_c*pulse_time*x_obs*R**.5)/(2*(math.pi*D*(t_obs**3))**.5))           *np.exp(-((R*x_obs-v*t_obs)**2)/(4*R*D*t_obs))*np.exp(-decay*t_obs) 


# In[12]:


plt.plot(t_obs,C[:,:])


# In[13]:


C=np.zeros((np.shape(t_obs)[0],np.shape(R)[0]))
  countR=-1
  for Rval in R:
      countR+=1
      C[:,countR]=((initial_c*pulse_time*x_obs*Rval**.5)/(2*(math.pi*D[countR]*(t_obs**3))**.5))*np.exp(-((Rval*x_obs-v[countR]*t_obs)**2)/(4*Rval*D[countR]*t_obs))*np.exp(-decay[countR]*t_obs)                          
  C[np.where(C<0)]=0
  C[np.where(np.isinf(C)==1)]=0
  C[np.where(np.isnan(C)==1)]=0
  for i in np.arange(trials):
    labelstr='run ' + str(i+1)
    plt.plot(t_obs,C[:,i],label=labelstr)
  plt.plot(t_obs,C_true,'.',label='C true w/error')
  plt.xlabel('Time (days)')
  plt.ylabel('Concentration')
  plt.title('X_obs =  200m') 
  #plt.legend()
  plt.show()


# In[14]:


numreparam=10
for j in np.arange(numreparam):
#for numreparam times, there will be nummetropolis parameter values generated in the same random uniform fashion as original
  nummetropolis = 5
    
    
    
  for k in np.arange(nummetropolis):
    v_metropolis=np.random.uniform((vmin),(vmax),size=1).round(4)
    D_metropolis=np.random.uniform((Dmin),(Dmax),size=1).round(4)
    R_metropolis=np.random.uniform((Rmin),(Rmax),size=1).round(4)
    decay_metropolis=np.random.uniform((decaymin),(decaymax),size=1).round(4)
    C_metropolis=np.zeros((np.shape(t_obs)[0]))
    C_metropolis[:]=((initial_c*pulse_time*x_obs*R_metropolis**.5)/(2*(math.pi*D_metropolis*(t_obs**3))**.5))*np.exp(-((R_metropolis*x_obs-v_metropolis*t_obs)**2)/(4*R_metropolis*D_metropolis*t_obs))*np.exp(-decay_metropolis*t_obs)                          
    C_metropolis[np.where(C_metropolis<0)]=0
    C_metropolis[np.where(np.isinf(C_metropolis)==1)]=0
    C_metropolis[np.where(np.isnan(C_metropolis)==1)]=0
    


    SE_metropolis = (C_metropolis-C_true)**2
    SSE_metropolis=np.sum(SE_metropolis)

    SSE=np.zeros((np.shape(R)[0]))
    countSSE=-1
    for i in np.arange(trials):
        SE = (C[:,i]-C_true)**2
        SSE[i]=np.sum(SE)      
    #within the numreparam times, runs with the highest SSE will be replaced with metropolis values
    worstrun=np.where(SSE==np.max(SSE))
    if (SSE[worstrun]>SSE_metropolis):# and (SSE[worstrun]>.5): #this leaves some original "good" values in
      v[worstrun]=v_metropolis 
      D[worstrun]=D_metropolis
      R[worstrun]=R_metropolis
      decay[worstrun]=decay_metropolis
      

                      
  L = 1/SSE
  L = L/np.sum(L)                          # likelihood weighted averaging
  
  Rnew=np.sum(R*L)
  Dnew=np.sum(D*L)
  decaynew=np.sum(decay*L)
  vnew=np.sum(v*L)
  
  worstrun=np.where(L==np.min(L))

  R[worstrun]=Rnew
  D[worstrun]=Dnew
  decay[worstrun]=decaynew
  v[worstrun]=vnew



  C=np.zeros((np.shape(t_obs)[0],np.shape(R)[0]))
  countR=-1
  for Rval in R:
      countR+=1
      C[:,countR]=((initial_c*pulse_time*x_obs*Rval**.5)/(2*(math.pi*D[countR]*(t_obs**3))**.5))*np.exp(-((Rval*x_obs-v[countR]*t_obs)**2)/(4*Rval*D[countR]*t_obs))*np.exp(-decay[countR]*t_obs)                          
  C[np.where(C<0)]=0
  C[np.where(np.isinf(C)==1)]=0
  C[np.where(np.isnan(C)==1)]=0
  for i in np.arange(trials):
    labelstr='run ' + str(i+1)
    plt.plot(t_obs,C[:,i],label=labelstr)
  plt.plot(t_obs,C_true,'.',label='C true w/error')
  plt.xlabel('Time (days)')
  plt.ylabel('Concentration')
  plt.title('X_obs =  200m') 
  #plt.legend()
  plt.show()


# In[15]:



trials=np.linspace(0,100,100)    #this loses our connection with the actual iterations

#maybe should be "for i in trials" with trials=100 at top of notebook
plt.plot(param_test1,trials,'.')
plt.xlabel('100 trials')
plt.ylabel('du and dD')
plt.title('Velocity and Dispersion space') 


# In[16]:


plt.plot(trials,param_test2,'.',color='green')
plt.xlabel('100 trials')
plt.ylabel('dR and ddecay')
plt.title('Retardation and Decay space') 


# In[ ]:





# In[17]:


plt.plot(u,u_orig,'.')
plt.xlabel('u final 100')
plt.ylabel('u original 100')
plt.title('u space') 
plt.show()

plt.plot(u_orig,'.')
plt.xlabel('u original 100')
plt.show()
plt.plot(umetro_sub,'.')
plt.xlabel('u metro sub')
plt.show()
plt.plot(u,'.')
plt.xlabel('u final 100')
plt.show()

x1=u_orig
x2=umetro_sub
x3=u
plt.plot(x1,'.')
#plt.plot(x2,'.')
plt.plot(x3,'.')
plt.xlabel('u original and u final')
plt.show()


# In[ ]:




