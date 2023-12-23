#!/usr/bin/env python
# coding: utf-8

#-----------Setting: Unsupervised fixed N, alpha=0.15------------

import numpy as np
import random
import pandas as pd
random.seed(10)

def which(self):
    try:
        self = list(iter(self))
    except TypeError as e:
        raise Exception("""'which' method can only be applied to iterables.
        {}""".format(str(e)))
    indices = [i for i, x in enumerate(self) if bool(x) == True]
    return(indices)


# In[30]:


K=20 #number of different groups


# In[31]:
m_len1=[]
m_ind1=[]
std_1=[]
m_len2=[]
m_ind2=[]
std_2=[]
m_len3=[]
m_ind3=[]
std_3=[]
m_len4=[]
m_ind4=[]
std_4=[]
m_len5=[]
m_ind5=[]
std_5=[]


for rep1 in range(40):

    ind1=[]
    ind2=[]
    ind3=[]
    ind4=[]
    ind5=[]
    len1=[]
    len2=[]
    len3=[]
    len4=[]
    len5=[]
    for rep in range(100):
        print([rep1,rep])
        s0=np.random.normal(0,10,K)
        N=[]
        for i in range(len(s0)):
            if s0[i]<0:
                N.append(15)  #generate 15 observations per tree
            else:
                N.append(15)
        s = np.random.normal(0,10,K) #mu
        mat=[]
        for j in range(K):
            mat.append(np.random.normal(s[j],0.5, N[j]))
        true=mat[K-1][-1]

        r=np.array(range(-1000,1000,1))/20
        validate=[0 for _ in range(len(r))]
        validate1=[0 for _ in range(len(r))]
        validate2=[0 for _ in range(len(r))]
        validate3=[0 for _ in range(len(r))]
        validate4=[0 for _ in range(len(r))]
        for i in range(len(r)):
            mat[K-1][-1]=r[i]
            mat_new=[]
            mean_val=[]
            mean_val2=[]
            for j in range(K):
                mean_val.append(np.mean(mat[j]))
                mean_val2.extend(mat[j])
            b_value=np.mean(mean_val)
            b_value1=np.mean(mean_val2)
            for j in range(K):
                #if N[j]==20:
                if np.sqrt(N[j])*abs(np.mean(mat[j])-b_value)/np.std(mat[j])<=2:
                    mat_new.extend(abs(mat[j]-b_value))
                else:
                    mat_new.extend(abs(mat[j]-np.mean(mat[j])))
            rk=sum(mat_new[-1]<mat_new)  #Our method
            if rk<45: 
                validate[i]=1 

            mat_new1=[]    #Lee's method  HCP (will only be used when the number of samples is random)
            for j in range(K):
                mat_new1.extend(abs(mat[j]-b_value))
            rk=sum(mat_new1[-1]<mat_new1)
            if rk<45:
                validate1[i]=1 

            mat_new2=[]     #conformal prediction method
            for j in range(K):
                mat_new2.extend(abs(mat[j]-b_value1))   
            rk=sum(mat_new2[-1]<mat_new2)
            if rk<round(sum(N)*0.15):
                validate2[i]=1
            
            mat_new3=[]     #sub-tree only
            mat_new3.extend(abs(mat[K-1]-np.mean(mat[K-1])))
            rk=sum(mat_new3[-1]<mat_new3)
            if rk<(np.floor(N[-1]*0.15)):
                validate3[i]=1
                
            mat_new41=[]     #Dunn's sub-sampling
            for j in range(K-1):
                mat_new41.append(mat[j][0])
            mat_new41.append(mat[K-1][-1])
            mat_new4=abs(np.array(mat_new41)-np.mean(mat_new41))
            rk=sum(mat_new4[-1]<mat_new4)
            if rk<(np.floor((K)*0.15)):
                validate4[i]=1
    #------------------Record the lengths of intervals-------------------
        len1.append(r[which(np.array(validate)==0)[-1]+1]-r[which(np.array(validate)==0)[0]])
        #np.savetxt("len1.csv", len1, delimiter=",")
        ind1.append(int(r[which(np.array(validate)==0)[0]]<=true<=r[which(np.array(validate)==0)[-1]+1]))
        #np.savetxt("ind1.csv", ind1, delimiter=",")
        len2.append(r[which(np.array(validate1)==0)[-1]+1]-r[which(np.array(validate1)==0)[0]])
        #np.savetxt("len2.csv", len2, delimiter=",")
        ind2.append(int(r[which(np.array(validate1)==0)[0]]<=true<=r[which(np.array(validate1)==0)[-1]+1]))
        #np.savetxt("ind2.csv", ind2, delimiter=",")
        len3.append(r[which(np.array(validate2)==0)[-1]+1]-r[which(np.array(validate2)==0)[0]])
        #np.savetxt("len3.csv", len3, delimiter=",")
        ind3.append(int(r[which(np.array(validate2)==0)[0]]<=true<=r[which(np.array(validate2)==0)[-1]+1]))
        #np.savetxt("ind3.csv", ind3, delimiter=",")
        len4.append(r[which(np.array(validate3)==0)[-1]+1]-r[which(np.array(validate3)==0)[0]])
        #np.savetxt("len4.csv", len4, delimiter=",")
        ind4.append(int(r[which(np.array(validate3)==0)[0]]<=true<=r[which(np.array(validate3)==0)[-1]+1]))
        #np.savetxt("ind4.csv", ind4, delimiter=",")
        len5.append(r[which(np.array(validate4)==0)[-1]+1]-r[which(np.array(validate4)==0)[0]])
        #np.savetxt("len5.csv", len5, delimiter=",")
        ind5.append(int(r[which(np.array(validate4)==0)[0]]<=true<=r[which(np.array(validate4)==0)[-1]+1]))
        #np.savetxt("ind5.csv", ind5, delimiter=",")

    m_len1.append(np.mean(len1))                     #Length of our method (SymmPI)
    np.savetxt("m_len1.csv",m_len1, delimiter=",")
    std_1.append(np.std(len1))                       #standard error of our method (SymmPI)
    np.savetxt("std1.csv",std_1, delimiter=",")
    m_ind1.append(np.mean(ind1))                     #coverage indicator of our method (SymmPI)
    np.savetxt("m_ind1.csv",m_ind1, delimiter=",")
    m_len2.append(np.mean(len2))                     #Length of HCP
    np.savetxt("m_len2.csv",m_len2, delimiter=",")
    std_2.append(np.std(len2))                       #standard error of HCP
    np.savetxt("std2.csv",std_2, delimiter=",")
    m_ind2.append(np.mean(ind2))                     #coverage indicator of HCP
    np.savetxt("m_ind2.csv",m_ind2, delimiter=",")
    m_len3.append(np.mean(len3))
    np.savetxt("m_len3.csv",m_len3, delimiter=",")    #Length of conformal
    std_3.append(np.std(len3))
    np.savetxt("std3.csv",std_3, delimiter=",")       #standard error of conformal
    m_ind3.append(np.mean(ind3))
    np.savetxt("m_ind3.csv",m_ind3, delimiter=",")    #coverage indicator of conformal
    m_len4.append(np.mean(len4))
    np.savetxt("m_len4.csv",m_len4, delimiter=",")   #length of sub-samples from last branch
    std_4.append(np.std(len4))
    np.savetxt("std4.csv",std_4, delimiter=",")      #standard error of sub-samples from last branch
    m_ind4.append(np.mean(ind4))
    np.savetxt("m_ind4.csv",m_ind4, delimiter=",")   #coverage indicator of sub-samples from last branch
    m_len5.append(np.mean(len5))                     #Length of sub-sampling
    np.savetxt("m_len5.csv",m_len5, delimiter=",")
    std_5.append(np.std(len5))
    np.savetxt("std5.csv",std_5, delimiter=",")     #standard error of sub-sampling
    m_ind5.append(np.mean(ind5))
    np.savetxt("m_ind5.csv",m_ind5, delimiter=",")   #coverage indicator of sub-sampling







