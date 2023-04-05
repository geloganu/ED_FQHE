import os

import csv
import numpy as np
import scipy.sparse
import scipy.linalg as LA

from sympy.physics.quantum.cg import CG

def nvec(orbs,Norb):
    ns = np.zeros(Norb)
    ns[orbs] = 1

    return ns

def I(ns):
    I = 0
    for i in range(0,len(ns)):
        I += (2**(i)*ns[i])
    return I
    
def m(alpha, ell):
    return alpha - ell

def Lplus(l, m):
    return np.sqrt(l*(l+1)-m*(m+1))

def Lminus(l, m):
    return np.sqrt(l*(l+1)-m*(m-1))

def LQN(matrix):
    return 1/2*(-1+np.sqrt(1+4*matrix))

def cg_coeffs(j1, j2):
    """
    args:
    j1 = j2 = J value of respective electrons' angular momenta
    """
    #triangular condition
    J = np.arange((j1 - j2), j1+j2 + 1)

    table = np.empty((int(2*j1)+1,2*int(j2+1),len(J)),dtype=object)

    for j1prime in range(0, int(2*j1)+1):
        for j2prime in range(0, int(2*j2)+1):
            for Jprime in range(0, len(J)):

                coeffs = np.array([])
                for M in range(-Jprime, Jprime+1):

                    val = float(CG(j1, m(j1prime, j1), j2, m(j2prime, j2), Jprime, M).doit())
                    coeffs = np.append(coeffs, val)

                table[j1prime, j2prime, Jprime] = coeffs

    return table

def L_spectrum(L2,hamil):
    energies, eigenstates=LA.eigh(hamil)
    L=np.array([])
    ene=np.array([])
    for j in range(0,len(energies)):
        angmomentum=np.matmul(np.matmul(eigenstates[:,j],L2),eigenstates[:,j].transpose())
        angmomentum=LQN(angmomentum)
        L=np.append(L,angmomentum)
        ene=np.append(ene,energies[j])

    ene=ene-min(ene)

    return L,ene

def overlap_Vectors(v1,v2):
    x=abs(np.dot(v1,v2)/np.linalg.norm(v1))
    print('overlap of wavefunciton <Ψtrial|Ψexact>=',x)
    return x

def get_eig(hamil):
    return LA.eigh(hamil.h.A)

def write_csv(array,file_name):
    with open(file_name,"w+") as fp:
        writer = csv.writer(fp,delimiter=",")
        writer.writerows(np.stack(array).T)

    