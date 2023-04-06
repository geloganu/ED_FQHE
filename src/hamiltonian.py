from itertools import combinations
from logging.config import valid_ident
import time
import numpy as np
from scipy.sparse import diags

from misc import *
from haldane_pseudopotential import *


class system:
    def __init__(self, N, Nphi, Lz, L2=False):
        """
        args:
        N: number of electrons
        Nphi: number of flux quantum
        Lz: angular momentum
        l:ell 
        Norb: number of orbitals
        Q: magnetic field stength?
        mzvals: list of m value
        occ_orbitals: list of occupied orbital states
        occ_orbitals_nlist: bit string of occupied orbitals
        occ_orbitals_ilist: integer list of occupid orbitals
        """

        #input constants
        self.N = N
        self.Nphi = Nphi
        self.Lz = Lz
        
        #defined constants
        self.Norb = Nphi + 1
        self.Q = Nphi/2
        self.mzvals = np.arange(-self.Q,self.Q+1)

        #orbital list
        occ_orbitals = np.zeros(self.N) 
        occ_orbitals_nlist = np.zeros(self.Norb)
        occ_orbitals_ilist = np.array([])

        #combinatorics all possible combination of states
        print('========Initializing possible states in spherical geometry========')
        st = time.time()
        for orbs in combinations(np.arange(0,self.Norb),self.N):
            orbs = np.array(orbs)

            #ensures states fulfill total angular momentum condition
            if sum(self.mzvals[orbs]) == Lz:
                occ_orbitals = np.vstack([occ_orbitals, orbs])
                
                ns = nvec(orbs,self.Norb)
                occ_orbitals_nlist = np.vstack([occ_orbitals_nlist, ns])
                
                occ_orbitals_ilist = np.append(occ_orbitals_ilist, I(ns))
        
        #delete initialization row in array
        occ_orbitals = np.delete(occ_orbitals,(0), axis = 0)
        occ_orbitals_nlist = np.delete(occ_orbitals_nlist,(0), axis = 0)

        #list of ind to sort lists
        ind = np.argsort(occ_orbitals_ilist)

        #sorted lists
        self.occ_orbitals = occ_orbitals[ind]
        self.occ_orbitals_nlist = occ_orbitals_nlist[ind]
        self.occ_orbitals_ilist = occ_orbitals_ilist[ind]
        self.sys_dim=len(self.occ_orbitals)
        
        print('completed in',time.time()-st,'seconds')
        print('N =',self.N,'electrons')
        print('Nphi =',self.Nphi,'flux quanta')
        print('')

        #initialize angular momentum matrix (need this?)
        if L2==True:
            self.L2 = self.angular_matrix()

    def get_ell(self):
        return self.Q

    def angular_matrix(self): 
        #L^2 = Lx^2 + Ly^2 + Lz^2 = 1/2(L+L- + L-L+) + Lz^2 using Lx^2 = 1/2(L+ + L-) + Ly^2 = (L+ - L-)
        
        #note indice change from alpha to x
        Q = self.Nphi/2
        M = self.Nphi + 1

        #L_plus L_minus constructor
        rows=np.array([])
        cols=np.array([])
        LpLmvals=np.array([])
        LmLpvals=np.array([])

        print('========Constructing total angular momentum matrix========')
        print(' ')
        st=time.time()

        for j in range(0,len(self.occ_orbitals)):
            ns = self.occ_orbitals_nlist[j]

            rows=np.append(rows,j)
            cols=np.append(cols,j)
            pmval=0
            mpval=0
            for alpha in range(0,M):
                pmval+=Lplus(Q,m(alpha,Q)-1)*Lminus(Q,m(alpha,Q))*ns[alpha]
                mpval+=Lminus(Q,m(alpha,Q)+1)*Lplus(Q,m(alpha,Q))*ns[alpha]
            LpLmvals=np.append(LpLmvals,pmval)
            LmLpvals=np.append(LmLpvals,mpval)

            for a1 in range(0,M-1):
                for a2 in range(1,M):
                    if a1!=a2 and (a1+1)!=(a2-1):
                        if ns[a1] == 1 and ns[a2] == 1:
                            sgn=1 if a1<a2 else -1
                            ns_new=ns.copy()
                            ns_new[a1]=0
                            ns_new[a2]=0

                            if ns_new[a1+1]==0 and ns_new[a2-1]==0:
                                sgn*=1 if (a1+1)<(a2-1) else -1
                                diagval=Lplus(Q,m(a1,Q))*Lminus(Q,m(a2,Q))

                                if diagval!=0:
                                    p1=sum(ns[0:a1])
                                    p2=sum(ns[0:a2]) 
                                    p1p1=sum(ns_new[0:(a1+1)])
                                    p2m1=sum(ns_new[0:(a2-1)])
                                    sgn*=(-1)**(p1 + p2 + p1p1 + p2m1)
                                    ns_new[a1+1]=1
                                    ns_new[a2-1]=1

                                    i = np.searchsorted(self.occ_orbitals_ilist, I(ns_new))

                                    rows=np.append(rows, i)
                                    cols=np.append(cols, j)
                                    LpLmvals=np.append(LpLmvals, -sgn*diagval)
                                    LmLpvals=np.append(LmLpvals, -sgn*diagval)

        LpLm = scipy.sparse.csc_matrix((LpLmvals, (rows, cols)), (self.sys_dim,self.sys_dim))
        LmLp = scipy.sparse.csc_matrix((LmLpvals, (rows, cols)), (self.sys_dim,self.sys_dim))

        #Lz^2 constructor
        Lz2vals = np.array([])
        for j in range(0,len(self.occ_orbitals)):
            val = 0
            for alpha in range(0,M):

                val += m(alpha,Q) * self.occ_orbitals_nlist[j][alpha]
                
            val = val*val
            Lz2vals = np.append(Lz2vals, val)

        Lz2 = diags(Lz2vals)
        
        L2=0.5*(LpLm+LmLp)+Lz2

        print('completed in',time.time()-st,'seconds')
        return L2
        
    def entanglement_spectrum(self,groundstate, subsystemA, NA, LzAvec):
        groundstate = groundstate/LA.norm(groundstate)
        dim=len(groundstate)

        IAlist = np.zeros(groundstate.shape[0], dtype=int)
        NAlist = np.zeros(groundstate.shape[0], dtype=int)
        LzAlist = np.zeros(groundstate.shape[0], dtype=float)
        mzA = self.mzvals[subsystemA]

        for i in range(0,dim):
            nsA = self.occ_orbitals_nlist[i][subsystemA]
            IAlist[i] = I(nsA)
            NAlist[i] = np.sum(nsA)
            LzAlist[i] = np.dot(mzA, nsA)
        
        print(f"Calculating entanglement spectrum for NA = {NA}, LzA = {LzAvec} ...")

        ent_spectrum = np.zeros((0, 2), dtype=float)

        for LzA in LzAvec:
            strip_inds = np.intersect1d(np.where(NAlist == NA), np.where(LzAlist == LzA))
            IAlist_stripped = IAlist[strip_inds]
            sort_inds = np.argsort(IAlist_stripped)
            IAlist_stripped = IAlist_stripped[sort_inds]

            groundstate_stripped = groundstate[strip_inds][sort_inds]

            nc = np.sum(IAlist_stripped == IAlist_stripped[0])
            nr = int(groundstate_stripped.shape[0] / nc)
            groundstate_matrix = groundstate_stripped.reshape((nr, nc)).T

            rhoA = np.dot(groundstate_matrix, groundstate_matrix.T)
            evals = LA.eigvals(rhoA)
            xi = -np.log(evals[np.where(evals > 0)])

            ent_spectrum = np.vstack((ent_spectrum, np.hstack((np.full((len(xi), 1), LzA), np.sort(xi).reshape(-1, 1)))))
        return ent_spectrum



        

class spherical_system:
    def __init__(self, system, pp_matrix):
        self.Norb = system.Norb
        self.nlist = system.occ_orbitals_nlist
        self.ilist = system.occ_orbitals_ilist
        self.sys_dim = len(system.occ_orbitals)
        self.pp_matrix = pp_matrix

        self.h = self.hconstructor()

    def hconstructor(self):
        dim = self.Norb

        rows = np.array([])
        cols = np.array([])
        hvals = np.array([])

        print('========Constructing Hamiltonian of size',self.sys_dim,'x',self.sys_dim,'========')
        st=time.time()

        for hcol in range(0,self.sys_dim):

            if hcol == 0 or hcol%np.ceil(0.05*self.sys_dim)==0:
                print("Working on interaction matrix:",hcol/self.sys_dim*100,'%')

            ns = self.nlist[hcol]
        
            for x2 in range(0, dim):
                for x1 in range(0, x2):
                    for x3 in range(0, dim):
                        for x4 in range(0, x3):
                            if x1 + x2 == x3+ x4:
                                if ns[x3] == 1 and ns[x4] == 1:
                                    ns_new = ns.copy()
                                    ns_new[x3] = 0
                                    ns_new[x4] = 0
                                    if ns_new[x1] == 0 and ns_new[x2] == 0:
                                        V1234 = self.pp_matrix[x1, x2, x3, x4]
                                        V1243 = self.pp_matrix[x1, x2, x4, x3]
                                        V2134 = self.pp_matrix[x2, x1, x3, x4]
                                        V2143 = self.pp_matrix[x2, x1, x4, x3]
                                        val = V1234 - V1243 - V2134 + V2143
                                        if val != 0:
                                            p1 = sum(ns_new[0:x1])
                                            p2 = sum(ns_new[0:x2])
                                            p3 = sum(ns[0:x3])
                                            p4 = sum(ns[0:x4])
                                            sgn = (-1)**(p1 + p2 + p3 + p4)
                                            ns_new[x1] = 1
                                            ns_new[x2] = 1
                                            i = np.searchsorted(self.ilist, I(ns_new))

                                            rows = np.append(rows, i)
                                            cols = np.append(cols, hcol)
                                            hvals = np.append(hvals, 0.5* sgn * val)

        hamiltonian = scipy.sparse.csc_matrix((hvals, (rows, cols)), (self.sys_dim,self.sys_dim))
        print('completed in',time.time()-st,'seconds')
        print('')
        return hamiltonian