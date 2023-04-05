import numpy as np
import matplotlib.pyplot as plt

def lowest_ene_spectrum(L,ene,lim):
    #ensure lowest energy eigenstates are selected

    ind=np.argsort(ene)
    sorted_L=L[ind]
    sorted_ene=ene[ind]
    if len(ene)>=lim: 
        sorted_L=sorted_L[:lim]
        sorted_ene=sorted_ene[:lim]
    
    return sorted_L, sorted_ene

def plot_L2_spectrum(L,ene,title,sorted=False,lim=125):
    #plot L2 spectrum with matplotlib
    if sorted==True:
        L,ene = lowest_ene_spectrum(L,ene,lim)

    ax=plt.plot(L, ene,ls="none", marker="_", ms="12", mew="1.5")
    plt.title(title)
    plt.xlabel('$L$')
    plt.ylabel('$E$')

    plt.show()
