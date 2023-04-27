import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as LA

from misc import *

def lowest_ene_spectrum(L,ene,lim):
    #ensure lowest energy eigenstates are selected

    ind=np.argsort(ene)
    sorted_L=L[ind]
    sorted_ene=ene[ind]
    if len(ene)>=lim: 
        sorted_L=sorted_L[:lim]
        sorted_ene=sorted_ene[:lim]
    
    return sorted_L, sorted_ene

def plot_L2_spectrum(hamil,L2,title,figsize,sorted=False,lim=125):
    #plot L2 spectrum with matplotlib
    L,ene=L_spectrum(L2,hamil)

    if sorted==True:
        L,ene = lowest_ene_spectrum(L,ene,lim)

    fig=plt.subplots(figsize=figsize)
    ax=plt.plot(L, ene,ls="none", marker="_", ms="12", mew="1.5")
    plt.title(title)
    plt.xlabel('$L$')
    plt.ylabel('$E$')

    plt.show()
