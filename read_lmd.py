#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 7 20:22:27 2022

@author: divya
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
import numpy as np
from statistics import mean
import random

# List initialization

MF1p,MF2p,MF1po,MF2po,KEF,Ef,chan = [],[],[],[],[],[],[]
MF1p_cor,MF2p_cor = [],[]

fig,ax=plt.subplots()

with open("data.lmd","r") as f:

# Parameter initialization    
    cha1=0;  cha2=0;  cha3=0;  cha4=0;  cha5=0;  cha6=0;
    ss1=0; ss2=0; ss3=0;  ss4=0; ss5=0;
    
    cha1_h=0;  cha2_h=0;  cha3_h=0;  cha4_h=0;  cha5_h=0;  cha6_h=0;  
    ss1_h=0; ss2_h=0; ss3_h=0;  ss4_h=0; ss5_h=0;
    
    cha1_s=0;  cha2_s=0;  cha3_s=0;  cha4_s=0;  cha5_s=0;  cha6_s=0;  
    ss1_s=0; ss2_s=0; ss3_s=0;  ss4_s=0; ss5_s=0;
    
    k=0;  j=0

    prob = {"ch1":"3000000000","ch2":"3300000000","ch3":"3330000000","ch4":"3333000000","ch5":"3333300000","ch6":"33333300"}
    sads = {"s1":"5000000000","s2":"5500000000","s3":"5550000000","s4":"5555000000","s5":"5555500000"}

    for line in f:                   
        mf1p = (line.split()[2])     
        mf2p = (line.split()[3])
        mf1po = (line.split()[4])
        mf2po = (line.split()[5])
        ke = (line.split()[21])
        particle     = (line.split())  

# Convoluting the values of  "mf1p" and "mf2p" events with the detection uncertainities    
        mf1p_val = 0.1*int(mf1p)
        val1=random.uniform(-mf1p_val,mf1p_val)
        mf1p_cor = int(mf1p)+val1
        
        mf2p_val = 0.045*int(mf2p)
        val2=random.uniform(-mf2p_val,mf2p_val)
        mf2p_cor = int(mf2p)+val2
        
        MF1p.append(mf1p)      
        MF2p.append(mf2p)
        MF1p_cor.append(mf1p_cor)      
        MF2p_cor.append(mf2p_cor)        
        MF1po.append(mf1po)
        MF2po.append(mf2po)
        KEF.append(ke)
        
        
#************  Total multiplicty w/o gate *************************

        for i in particle:
                 if (i==prob["ch1"]):cha1+=1
                 if (i==prob["ch2"]):cha2+=1
                 if(i==prob["ch3"]):cha3+=1
                 if(i==prob["ch4"]):cha4+=1
                 if(i==prob["ch5"]):cha5+=1
                 if(i==prob["ch6"]):cha6+=1

                 if(i==sads["s1"]): ss1+=1  
                 if(i==sads["s2"]):ss2+=1
                 if(i==sads["s3"]):ss3+=1
                 if(i==sads["s4"]):ss4+=1
                 if(i==sads["s5"]):ss5+=1
                 
#************  Gated condition on "H"   *************************
      
# Extract events with gated condition on "mf2p_cor"       
 
        if( int(mf2p_cor)>130 and int(mf2p_cor) < 144):
            k+=1
            for i in particle:
                
                if (i==prob["ch1"]):cha1_h+=1
                if(i==prob["ch2"]):cha2_h+=1
                if(i==prob["ch3"]):cha3_h+=1
                if(i==prob["ch4"]):cha4_h+=1
                if(i==prob["ch5"]):cha5_h+=1
                if(i==prob["ch6"]):cha6_h+=1

                if(i==sads["s1"]): ss1_h+=1  
                if(i==sads["s2"]):ss2_h+=1
                if(i==sads["s3"]):ss3_h+=1
                if(i==sads["s4"]):ss4_h+=1
                if(i==sads["s5"]):ss5_h+=1  

#************  Gated condition on "S"  *************************
 
# Extract events with gated condition on "mf1p_cor" and "mf2p_cor"        
                
        if( int(mf1p_cor)>=117 and int(mf2p_cor) <=130):
            
            j+=1
            for i in particle:
                
                if (i==prob["ch1"]):cha1_s+=1
                if(i==prob["ch2"]):cha2_s+=1
                if(i==prob["ch3"]):cha3_s+=1
                if(i==prob["ch4"]):cha4_s+=1
                if(i==prob["ch5"]):cha5_s+=1
                if(i==prob["ch6"]):cha6_s+=1

                if(i==sads["s1"]):ss1_s+=1  
                if(i==sads["s2"]):ss2_s+=1
                if(i==sads["s3"]):ss3_s+=1
                if(i==sads["s4"]):ss4_s+=1
                if(i==sads["s5"]):ss5_s+=1                           
 
#************  Distribution function w/o and with inclusion  of uncertainity in the data *************************

TM=MF1p+MF2p 
TM_cor=MF1p_cor+MF2p_cor    

ax=sn.distplot(TM,87,label="w/o uncertainity")   
ax=sn.distplot(TM_cor,89,label="w uncertainity")

ax.legend(loc="center",bbox_to_anchor=(0.8,0.85),fontsize=10)
ax.set_xlabel(r"FM",fontsize=10)
ax.set_ylabel(r"Normalized counts",fontsize=10)

#************  Calculating multiplicty value for entity "H" *************************

np_Hp=((cha1_h+2*cha2_h+3*cha3_h+4*cha4_h+5*cha5_h+6*cha6_h)/k)
np_Hpo=((ss1_h+2*ss2_h+3*ss3_h+4*ss4_h+5*ss5_h)/k)

print("total value of entity H :",np_Ap+np_Apo)

#************  Calculating multiplicty value for entity "S"  ************************


np_Sp=((cha1_s+2*cha2_s+3*cha3_s+4*cha4_s+5*cha5_s+6*cha6_s)/j)
np_Spo=((ss1_s+2*ss2_s+3*ss3_s+4*ss4_s+5*ss5_s)/j)

print("total value of entity S :",np_Sp+np_Spo)

#************ Calculating the total multiplicty value w/o gate  *********************

np_p=(chan1+2*chan2+3*chan3+4*chan4+5*chan5+6*chan6)/100000
np_po=(neuss1+2*neuss2+3*neuss3+4*neuss4+5*neuss5)/100000  

print("Total value w/o gate condition :",np_p+np_po)

#-------------------------------- END OF SCRIPT -------------------------------------