import time
import math
import numpy as np
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j

from misc import *
from hamiltonian import *

class haldane_pseudopotential:
    def __init__(self, l, LLn, custom=None):
        """
        args:
        l: angular momentum in the effective LLL problem == monopole strength Q in LLL
        LLn: nth Landau Level
        trial= (opt) trial pseudopotential
        """
        
        self.l = l
        self.LLn = LLn

        self.L = (np.arange(0, int(2*self.l) + 1,dtype=float)) #+1 for inclusive endpoint
        self.m = (2*self.l - self.L)
        
        ind=np.argsort(self.m)
        self.L=self.L[ind]
        self.m=self.m[ind]
        
        if custom is None:
            self.pseudopotential()
            self.pp_matrix = self.pp_matrix_generator(self.V)
        elif custom is not None:
            if isinstance(custom,np.ndarray): custom=np.asarray(custom)
            self.V_trial=custom
            self.pp_matrix=self.pp_matrix_generator(self.V_trial)


    def pseudopotential(self):
        """
        consts:
        L: total angular momentum on the sphere
        m: 2*l - L = relative momenta on sphere
        V: pseudopotential values
        """

        #local attribute for performance
        l=self.l 
        L=self.L

        Q = l - self.LLn
        
        print('========Initializing two-body pseudopotential========')
        st=time.time()

        self.V = np.zeros(len(L))

        Vk = np.full(int(2*l)+1,1/np.sqrt(l)) #radius=sqrt(self.l)
        

        #two body interaction value
        for i in range(0,len(L)):
            vk = 0

            #summation of Vk over k
            for k in range(0, int(2*l) + 1):
                vk += Vk[k] * wigner_6j(L[i], l, l, k, l, l) * wigner_3j(l, k ,l, -Q, 0, Q)**2
            
            self.V[i] = vk * (-1)**(2*Q + L[i]) * (2*l + 1)**2
        

        self.V = self.V[np.argsort(self.m)]

        print('completed in',time.time()-st,'seconds')
        print('pseudopotential:',self.V)
        print('')


    def pp_matrix_generator(self,pp):
        
        print('========Initializing interaction matrix========')
        cg_table = cg_coeffs(self.l, self.l)

        Norb = len(pp)
        pp_matrix = np.empty((Norb, Norb, Norb, Norb))
        st=time.time()
        for x1 in range(0,Norb):

            if x1 == 0 or x1%np.ceil(0.05*Norb)==0:
                print("Working on interaction matrix:",x1/Norb*100,'%')
            
            for x2 in range(0,Norb):
                for x3 in range(0,Norb):
                    for x4 in range(0,Norb):
                        V = 0

                        for l in range(0, int(2*self.l)+1):
                            for Mz in range(-int(l),int(l+1)):
                                V += cg_table[x1, x2, l][int(Mz + l)]*pp[int(-1-l)]*np.conj(cg_table[x3, x4, l][int(Mz + l)])

                                
                        
                        pp_matrix[x1, x2, x3, x4] = V
        
        print('completed in ',time.time()-st,'seconds')      
        print('')
        
        return pp_matrix

    def overview(self):
        print("---------------------")
        print("total angular momentum:", self.L)
        print('relative momenta:', self.m)
        print('pseudopotentials:', self.V)
        print("---------------------")