import numpy as np
import scipy.linalg as LA
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, '../src')

from hamiltonian import *
from misc import *
from haldane_pseudopotential import *

#simulation parameters
#for MR_Pf state, v=5/2 and S=3
#N_phi=2(N_e)-3 for LL=2
N = 10
filling_factor_inv = 2
top_shift=3
Nphi = filling_factor_inv*N-top_shift


#initialization (system configurations and parameters)
system = system(N, Nphi, 0)
#print(system.Q)
#print(system.L2.A)

#trial pseudopotentials
trial_pp=np.array([0,1.4,0,1])
trial_pp=np.append(trial_pp, np.zeros((Nphi+1)-len(trial_pp)))

#Coulomb pseudopotentials
pp = haldane_pseudopotential(system.Q, LLn = 1)

#trial pseudopotentials
pp_trial = haldane_pseudopotential(system.Q, LLn = 1,custom=trial_pp)
#pp1 = haldane_pseudopotential(system.Q, LLn = 1)
#pp.overview()

#create hamiltonian and diagonalize coulomb interaction
hamiltonian = spherical_system(system, pp.pp_matrix)
energies, eigenstates=get_eig(hamiltonian)

#calculate trial wavefunction using modified Haldane pseudopotential terms
trial_hamiltonian = spherical_system(system, pp_trial.pp_matrix)
trial_energies, trial_eigenstates=get_eig(trial_hamiltonian)

#overlap of wavefunction output
overlap_Vectors(np.transpose(trial_eigenstates)[0],np.transpose(eigenstates)[0])

