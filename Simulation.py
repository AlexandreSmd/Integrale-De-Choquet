# -*- coding: utf-8 -*-
"""
Created on Sun May 15 16:41:33 2022

@author: asima
"""



import numpy as np
import Capacity
import random as rd
from math import *
import matplotlib.pyplot as plt
from numpy import log

from scipy.stats import beta, gaussian_kde

c=2
cc=5

def powerset(original_set):
  # below gives us a set with one empty set in it
  ps = set({frozenset()}) 
  for member in original_set:
    subset = set()
    for m in ps:
      # to be added into subset, needs to be
      # frozenset.union(set) so it's hashable
      subset.add(m.union(set([member])))
    ps = ps.union(subset)
  return ps


# Cas avec deux aléa

def capaciteMu2(n):
    capa=powerset(set([i for i in range(1,n+1)])) #on a tous les sous-ensembles 
    #print(capa)
    keys=[]
    values=[]
    for key in capa:    #on créer a liste des sous ensembles
        keys.append(key)
    keys.sort(key=len)      #on a ainsi la liste des sous-ensembles trier par ordre croissants de leur dimension
    
    #m=len(keys)
    capacite={keys[0]:0}    #on créer un dictionnaire
    values=[0,0.07,0.2,0.1,0.15,0.5,0.5,0.6,0.35,0.44,0.4,0.8,0.9,0.9,0.76,1]
    C1=rd.uniform(0.2,0.8)
    C2=C1 + np.random.beta(c,cc)
    if C2>=1 or C2<=values[6] or C2 <=values[10]:
        C2=(max(C1,values[6],values[10])+1)/2
    values[5]=C1
    values[14]=C2
    a=C1
    b=C2
    
    for i in range(1,16):
        capacite[keys[i]]=values[i]
    
    capaciteF = Capacity.Capacity(capacite)
    print(capacite)
    return [capaciteF,a,b]



################################"
##Affichage grapiques des résultats : 
    
def graphe():
    M=0   #liste des intégrales de Choquet avec un coefficient aléatoire
#exemple de constante de capacité
    X=[2,1,7,4]
    X=np.reshape(X, (1, 4))
    L=[]
    coeff1=[]
    coeff2=[]

    somme=0
    somme2=0
    for i in range(100000):
        tirage=capaciteMu2(4)            
        capa=tirage[0]
        coeff1.append(tirage[1])
        coeff2.append(tirage[2])
        M=capa.CIntegral(X)
        L.append(M[0])
        somme+=M[0]
        somme2+=M[0]**2
    esp=somme/100000
    var=somme2/100000-esp**2
    ect=(var/100000)**0.5
    print(esp,'+/-',2*ect)
    print('X=',X)    
    return(L)
    
#########################################"
##FONCTION DE REPARTITION : 
    
def fonction_repartition(L):
    X=[2,1,7,4]
    X=np.reshape(X, (1, 4))
    n=len(L)
    x=np.linspace(2.5,4.7,100)
    F=[]
    f=0
    
    for i in range(100-1):
        
        for j in range(n):
            if L[j]>=x[i] and L[j]<=x[i+1]:
                f+=1
            
        F.append(f/n)
    F.append(1)
    plt.plot(x,F)
    plt.hist(L,100,cumulative=True,density=True)
    plt.title('fonction de répartition')
    plt.xlabel("x")
    plt.ylabel("F(Choquet)")
    
def repartition():
    X=[2,1,7,4]
    X=np.reshape(X, (1, 4))
    n=1000 #Nombre d'itérations générale et donc pas de notre fonction
    m=2000 #Raffinage du calcul de la proba
    x=np.linspace(2.7,4.5,n)
    F=[]
    f=0
    for i in range (n):
        somme=0
        
        for j in range(m):
            l=0
            tirage=capaciteMu2(4)            
            capa=tirage[0]
            M=capa.CIntegral(X)
            l=M[0]
            
            if float(l)<x[i]:
                somme+=1
        proba=somme/m
        f=proba
        F.append(f)
    return(x,F)

    



###MEILLEUR ALTERNATIVE POUR LE CALCUL DE LA DENSITÉ###

def density_scipy():
    data = graphe()
    density = gaussian_kde(data)
    x = np.linspace(2.5,5,300)
    y=density(x)

    plt.plot(x, y)
    plt.title("Density Plot of the data")
    plt.show()
