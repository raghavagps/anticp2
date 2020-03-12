import os
import numpy as np
def aac_gen(seq,option,x,y):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    seq = seq.upper()
    aac=[]
    if option=='Normal':
        seq=seq
    elif option=='N':
        seq=seq[0:x]
    elif option=='C':
        seq=seq[-x:][::-1]
    elif option=='NC':
        seq=seq[0:x]+seq[-y:][::-1]
    for i in std:
        counter = seq.count(i) 
        aac+=[((counter*1.0)/len(seq))*100]
    return aac            
        
def dpc_gen(seq,option,x,y):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    seq=seq.upper()
    dpc=[]
    if option=='Normal':
        seq=seq
    elif option=='N':
        seq=seq[0:x]
    elif option=='C':
        seq=seq[-x:0][::-1]
    elif option=='NC':
        seq=seq[0:x]+seq[-y:][::-1]
    for j in std:
        for k in std:
            temp  = j+k
            count = seq.count(temp)
            dpc+=[((count*1.0)/(len(seq)-1))*100]
    return dpc
            
def bin_aac(seq,option,x=None,y=None):
    amino_acids=list("ACDEFGHIKLMNPQRSTVWY")
    Dict={}
    #print(x)
    for i,j in enumerate(amino_acids):
        Dict[j]=i
    seq=seq.upper()
    lis=np.asarray([])
    if option=='Normal':
        seq=seq
    elif option=='N':
        seq=seq[0:x]
    elif option=='C':
        seq=seq[len(seq)-x:][::-1]
    elif option=='NC':
        seq=seq[0:x]+seq[len(seq)-y:][::-1]
    for i in seq:
        a=np.zeros(len(amino_acids))    
        a[Dict[i]]=1
        lis=np.append(lis,a)
    return lis
