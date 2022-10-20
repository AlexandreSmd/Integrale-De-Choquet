# -*- coding: utf-8 -*-
"""
Created on Tue May 10 11:21:59 2022

@author: asima
"""

import itertools
import operator
import numpy as np


class Capacity:
#Class representing a capacity, i.e. a monotone set function vanishing at the empty set (also called
#fuzzy measure, non-additive measure, monotone measure).
 
    def __init__(self, capacity_set):  #création d'un objet Capacity

        # vérification des propriété de la capacité:
        # 1- nombre d'éléments
        if np.log2(len(capacity_set)).astype(int) != np.log2(len(capacity_set)):
            print(f'The number of element in the Capacity is not consistent, i.e. 2^N')
            return 

 #Cardinalité de l'ensemble, S, sur le quel est définie la capacité
        self.N = np.log2(len(capacity_set)).astype(int)

        
#powerset de la capacité
        Set_S = list(range(1,self.N+1))
        PS = [[]]
        for x in Set_S:
            # for every additional element in our set the power set consists of the subsets that don't
            # contain this element (just take the previous power set) plus the subsets that do contain 
            # the element (use list comprehension to add [x] onto everything in the previous power set)
            PS.extend([subset + [x] for subset in PS])

        self.powerset = list(frozenset(i) for i in PS)


        # 2-Borne de la capacité
        if capacity_set[frozenset({})]!=0:
            print(f'The capacity of the empty set is not equal to 0.')
            return

        Set_S = list(range(1,self.N+1))
        if capacity_set[frozenset(Set_S)]!=1 :
            print(f'The capacity of the {Set_S} set is not equal to 1.')
            return


        # 3-Monotonicity
        for SS in self.powerset:
            for E in set(Set_S).difference(set(SS)):
                if (capacity_set[SS] <= capacity_set[SS.union(frozenset({E}))]) != True:
                    print(f'The capacity is not monotone.')
                    print(f'mu{SS} > mu{SS.union(frozenset({E}))} ---> {capacity_set[SS]} > {capacity_set[SS.union(frozenset({E}))]}')
                    return

       
        #Capacité elle même
        self.mu = capacity_set
        
#=================================================================================================
# AFFICHAGE DE LA CAPACITE
    def printC(self,Order="Binary"):
        #affiche les parties de S avec la valeur de la capacité associée
        print(f'================================================================')
        if Order == "Binary":
            for s in self.powerset:
                print(f'{set(s)} --> {self.mu[s]}')
        elif Order == "Natural":
            for s in sorted(self.powerset, key=len):
                print(f'{set(s)} --> {self.mu[s]}')
        else:
            print("Valeur de Type inconnue")
        print(f'================================================================')
    

#=================================================================================================
# CALCUL DE L'INTEGRALE DE CHOQUET

    def CIntegral(self,X):
        '''Calculate the Choquet Integral and return just that value'''
        '''if X is a vector, a scalar value is return, if X is a matrix, a vector is return'''

        X_sorted = np.sort(X)
        X_sorted_idx = np.argsort(X)
        X_sorted_diff=X_sorted[:,0]


        for i in range(1,self.N): # on calcul X(i)-X(i-1) dans une matrice
            X_sorted_diff= np.c_[X_sorted_diff,X_sorted[:,i]-X_sorted[:,i-1]]
 
        #initialisation de l'integrale
        CI = np.zeros(X.shape[0])
        for line in range(0,X.shape[0]): # on parcourt toutes les lignes
            for j in range(0,self.N): # on parcours chaque (i)
                CI[line]= CI[line] + X_sorted_diff[line,j]*self.mu[frozenset(X_sorted_idx[line,j:]+1)]
 
        return CI

#essais capacity



import numpy as np

import numpy as np
import random as rd
from math import *
import matplotlib.pyplot as plt
from numpy import log

from scipy.stats import beta, gaussian_kde

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

def exp(lam):
    u = rd.random()
    return - (1/lam) * log (1-u)


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
    # C2 = C1 + exp(8) #On définit nos variables aléatoires pour les coefs       
    # while C1 < 0.2 or C1 > 0.8:
    #     C1=exp(0.75)
    C2 = rd.uniform(C1, 1)
    while C2 < C1 or C2 > 1 or C2 < 0.65 or C2 < 0.4:
          C2=rd.uniform(C1, 1)
    #print(C1,C2)
    
    # while C2<0.7: #On retire tant que C2 n'est pas supérieur à 0.7
    #     C2=(float(beta.rvs(26, 2, loc=0, scale=1, size=1, random_state=None))-0.7)/0.3
    # while C1>0.5: #On retire tant que C1 n'est pas inférieur à 0.5
    #     C1=float(beta.rvs(3, 18, loc=0, scale=1, size=1, random_state=None))/0.5
    
    values[5]=C1
    values[14]=C2
    
    # values[14]=rd.uniform(0.7, 1)
    # values[2]=rd.uniform(0, 0.5)
    
    a=values[14]
    b=values[5]
    #for i in range(1,n):  #on parcours les sous ensembles suivant leur dimensions
        # c=int(fact(n)/(fact(i)*fact(n-i)))   #nb de parties à i elements dans un ensemble de n elements
        # c_1=int(fact(n)/(fact(i-1)*fact(n-i+1))) #nb de parties à i-1 elements dans un ensemble de n elements
        # l=int(sommeCombi(n,i-1)) #total d'ensemble de dimensions inférieure à i-1
        # for q in range(l,l+c):  #on itère seulement sur les ensemble de dimension i, donc ceux entre l et c+l
        #     SubsetsV=[]     #on créer une liste des sous-ensembles de keys[q]
        #     for j in range (l-c_1,l):     #on parcours les sous-ensembles de dimension i-1
        #         if keys[j].issubset(keys[q])==True:
        #             SubsetsV.append(values[j])
        #     M=maximum(SubsetsV)
        #     #v=np.random.uniform(SubsetsV[M],1)  #on vérifie ainsi la propriete de croissance de la capacité
        #     v=SubsetsV[M]+0.007*(q-l+1)        #pour n=7
        #     #v=SubsetsV[M]+0.02*(q-l+1)         #pour n=4
        #     if i==n-2 and q==l+c-1:
        #         v=np.random.uniform(SubsetsV[M],1)      #pour n=4 c'est le coeff {1,3} pour n=7 c'est {2,3,5,6,7}
        #     if q==m-2:
        #         v=np.random.uniform(SubsetsV[M],1)  #pour n=7 c'est le coeff {2,3,4,5,6,7} et pour n=4 c'est {1,3,4}
        #     values.append(v)
        #     capacite[keys[q]]=values[q] #on ajoute au dictionnaire le nouveau sous ensemble et sa valeur associée
        
    # capacite[keys[m-1]]=1 #on termine avec l'ensemble lui même 
    # values.append(1)
    for i in range(1,16):
        capacite[keys[i]]=values[i]
    
    capaciteF = Capacity.Capacity(capacite)
    print(capacite)
    return [capaciteF,b,a]



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
    for i in range(10000):
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
    
def grapheL(A):
    plt.hist(A,edgecolor='red',density=True)
    plt.title("Intégrale Choquet(X)")
def graphecoeff1(A):
    plt.hist(A,edgecolor='red',density=True)
    plt.title("mu({3})")
def graphecoeff2(A):
    plt.hist(A,edgecolor='red',density=True)
    plt.title("mu({1,3,4})")
    
#########################################"
##FONCTION DE REPARTITION : 
    
def fonction_repartition():
    X=[2,1,7,4]
    X=np.reshape(X, (1, 4))
    n=1000 #Nombre d'itérations générale et donc pas de notre fonction
    m=1000 #Raffinage du calcul de la proba
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
    plt.plot(x,F)
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

    

def densite(x,F):   # x = liste des abcisses, F = fonction de répartition
    n = 100
    L = []
    xnew = np.zeros(n-1)
    for k in range (0,n-1):
        xnew[k] = (x[10*k] + x[10*(k+1)]) / 2
        accroissement = (F[10*(k+1)] - F[10*k]) / (x[10*(k+1)] - x[10*k])   #calcul du taux d'accroissement
        L.append(accroissement)
    #L.append(L[-1])
    plt.plot(xnew,L)
    plt.title('densité')
    plt.xlabel("x")
    plt.ylabel("Fonction densité")

def densite2(x,F):   # x = liste des abcisses, F = fonction de répartition
    n = 10  #oit être un diviseur du nombre de points de la fct de répartition
    alpha = int(len(x)/n)
    print(alpha)
    L = []
    xnew = np.zeros(n-1)
    for k in range (0,n-1):
        xnew[k] = (x[alpha*k] + x[alpha*(k+1)]) / 2
        accroissement = (F[alpha*(k+1)] - F[alpha*k]) / (x[alpha*(k+1)] - x[alpha*k])   #calcul du taux d'accroissement
        L.append(accroissement)
    #L.append(L[-1])
    plt.plot(xnew,L)
    plt.title('densité')
    plt.xlabel("x")
    plt.ylabel("Fonction densité")
    return(xnew,L)


def densite_moy(x,F):   # x = liste des abcisses, F = fonction de répartition
    n = len(x)
    L = []
    xnew = np.zeros(n-1)
    for k in range (0,n-1):
        xnew[k] = (x[k] + x[k+1]) / 2
        if k <10 or k>n-11:
            accroissement = (F[k+1] - F[k]) / (x[k+1] - x[k])   #calcul du taux d'accroissement
            L.append(accroissement)
        else : 
            L_inter=[]
            for i in range(19):
                accroissement_inter = (F[k-10+i+1] - F[k-10+i]) / (x[k-10+i+1] - x[k-10+i])   #calcul du taux d'accroissement
                if accroissement_inter < 0:
                    L_inter.append(0)
                else :
                    L_inter.append(accroissement_inter)
            accroissement = np.mean(L_inter)
            L.append(accroissement)
    #L.append(L[-1])
    plt.plot(xnew,L)
    plt.title('densité')
    plt.xlabel("x")
    plt.ylabel("Fonction densité")

def rectangle(xnew,L):
    S = 0
    for k in range (len(L)-1):
        S+= ( (L[k+1] + L[k]) / 2 ) * (xnew[k+1] - xnew[k])
    return(S)


###MEILLEUR ALTERNATIVE POUR LE CALCUL DE LA DENSITÉ###

def density_scipy():
    data = graphe()
    density = gaussian_kde(data)
    x = np.linspace(2.5,5,300)
    y=density(x)

    plt.plot(x, y)
    plt.title("Density Plot of the data")
    plt.show()
