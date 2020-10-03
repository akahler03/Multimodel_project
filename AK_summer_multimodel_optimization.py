#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
initial_c=1
x_obs=200
pulse_time=200
t_start=0               
t_end = 800
dt = 10
t_obs=np.arange(t_start,t_end,dt)  


# In[2]:


trials=100   #np.linspace(0,100,100)
umin=1
umax=3
Dmin=10
Dmax=100
Rmin=1
Rmax=2
decaymin=0
decaymax=.005


# In[3]:


u=np.random.uniform((umin),(umax),size=trials).round(4)
D=np.random.uniform((Dmin),(Dmax),size=trials).round(4)
R=np.random.uniform((Rmin),(Rmax),size=trials).round(4)
decay=np.random.uniform((decaymin),(decaymax),size=trials).round(4)

utrue=1.6
Dtrue=55
Rtrue=1.3
decaytrue=0


# In[4]:


u_orig=u.copy()
D_orig=D.copy()
R_orig=R.copy()
decay_orig=decay.copy()
param_orig=[u_orig,D_orig,R_orig,decay_orig]


# In[5]:


C=np.zeros((np.shape(t_obs)[0],np.shape(R)[0]))
countR=-1
for Rval in R:
    countR+=1
    C[:,countR]=((initial_c*pulse_time*x_obs*Rval**.5)/(2*(math.pi*D[countR]*(t_obs**3))**.5))*np.exp(-((Rval*x_obs-u[countR]*t_obs)**2)/(4*Rval*D[countR]*t_obs))*np.exp(-decay[countR]*t_obs)                          
C[np.where(C<0)]=0
C[np.where(np.isinf(C)==1)]=0
C[np.where(np.isnan(C)==1)]=0

C_true=np.zeros((np.shape(t_obs)[0]))
C_true[:]=((initial_c*pulse_time*x_obs*Rtrue**.5)/(2*(math.pi*Dtrue*(t_obs**3))**.5))*np.exp(-((Rtrue*x_obs-utrue*t_obs)**2)/(4*Rtrue*Dtrue*t_obs))*np.exp(-decaytrue*t_obs)                          
C_true[np.where(C_true<0)]=0
C_true[np.where(np.isinf(C_true)==1)]=0
C_true[np.where(np.isnan(C_true)==1)]=0

err_val=np.random.normal(loc=0,scale=.03,size=(C_true.size))
C_true=C_true+err_val


# In[6]:


for i in np.arange(trials):
  labelstr='run ' + str(i+1)
  plt.plot(t_obs,C[:,i],label=labelstr)
plt.plot(t_obs,C_true,'.',label='C true w/error')
plt.xlabel('Time (days)')
plt.ylabel('Concentration')
plt.title('X_obs =  200m') 
#plt.legend()
plt.show()


# In[7]:


#keeps track of the chosen metropolis
umetro_sub=[]
Dmetro_sub=[]
Rmetro_sub=[]
decaymetro_sub=[]


numreparam = 80
for j in np.arange(numreparam):
#for numreparam times, there will be nummetropolis parameter values generated in the same random uniform fashion as original
  nummetropolis = 50
  for k in np.arange(nummetropolis):
    u_metropolis=np.random.uniform((umin),(umax),size=1).round(4)
    D_metropolis=np.random.uniform((Dmin),(Dmax),size=1).round(4)
    R_metropolis=np.random.uniform((Rmin),(Rmax),size=1).round(4)
    decay_metropolis=np.random.uniform((decaymin),(decaymax),size=1).round(4)
    C_metropolis=np.zeros((np.shape(t_obs)[0]))
    C_metropolis[:]=((initial_c*pulse_time*x_obs*R_metropolis**.5)/(2*(math.pi*D_metropolis*(t_obs**3))**.5))*np.exp(-((R_metropolis*x_obs-u_metropolis*t_obs)**2)/(4*R_metropolis*D_metropolis*t_obs))*np.exp(-decay_metropolis*t_obs)                          
    C_metropolis[np.where(C_metropolis<0)]=0
    C_metropolis[np.where(np.isinf(C_metropolis)==1)]=0
    C_metropolis[np.where(np.isnan(C_metropolis)==1)]=0
    
#NEW STUFF BELOW

#framework for distance calculation
    du=(u_metropolis-u_orig)**2
    dD=(D_metropolis-D_orig)**2
    dR=(R_metropolis-R_orig)**2
    ddecay=(decay_metropolis-decay_orig)**2
    space_rep=(du+dD+dR+ddecay)**1/2
   #Goal: force du to increase--to look further
   # test_space=5
    #for t in test_space:
     # something to adjust the parameters based on the value of space_rep 
        
    
    #NEW STUFF ABOVE
    
    
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
      u[worstrun]=u_metropolis 
      D[worstrun]=D_metropolis
      R[worstrun]=R_metropolis
      decay[worstrun]=decay_metropolis
      umetro_sub.append(u_metropolis[0])  
      Dmetro_sub.append(u_metropolis[0])
      Rmetro_sub.append(R_metropolis[0])
      decaymetro_sub.append(decay_metropolis[0])
    
   
        
  
    

                      
  L = 1/SSE
  L = L/np.sum(L)                          # likelihood weighted averaging
  
  Rnew=np.sum(R*L)
  Dnew=np.sum(D*L)
  decaynew=np.sum(decay*L)
  unew=np.sum(u*L)
  
  worstrun=np.where(L==np.min(L))

  R[worstrun]=Rnew
  D[worstrun]=Dnew
  decay[worstrun]=decaynew
  u[worstrun]=unew

  #NEW STUFF BELOW
  
    
    #NEW STUFF ABOVE

  C=np.zeros((np.shape(t_obs)[0],np.shape(R)[0]))
  countR=-1
  for Rval in R:
      countR+=1
      C[:,countR]=((initial_c*pulse_time*x_obs*Rval**.5)/(2*(math.pi*D[countR]*(t_obs**3))**.5))*np.exp(-((Rval*x_obs-u[countR]*t_obs)**2)/(4*Rval*D[countR]*t_obs))*np.exp(-decay[countR]*t_obs)                          
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


# In[8]:


min(du)


# In[9]:


max(du)


# In[8]:


du_space=max(du)-min(du)


# In[ ]:


du_space


# In[ ]:


print(umin)
print(umax)
print(umax-umin)


# In[9]:


du=(u_metropolis-u_orig)**2            #Condsidering ALL metropolis generated and the original 100
                                       #or should it only be CHOSEN metropolis?
dD=(D_metropolis-D_orig)**2
dR=(R_metropolis-R_orig)**2
ddecay=(decay_metropolis-decay_orig)**2
space_rep=(du+dD+dR+ddecay)**1/2
print(np.shape(space_rep))


# In[10]:


param_test1=(du+dD)**1/2
param_test2=(dR+ddecay)**1/2


# In[11]:



trials=np.linspace(0,100,100)    #this loses our connection with the actual iterations

#maybe should be "for i in trials" with trials=100 at top of notebook
plt.plot(param_test1,trials,'.')
plt.xlabel('100 trials')
plt.ylabel('du and dD')
plt.title('Velocity and Dispersion space') 


# In[11]:


plt.plot(trials,param_test2,'.',color='green')
plt.xlabel('100 trials')
plt.ylabel('dR and ddecay')
plt.title('Retardation and Decay space') 


# In[ ]:





# In[ ]:


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




