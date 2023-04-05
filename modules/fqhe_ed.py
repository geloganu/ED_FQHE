#!/usr/bin/env python

import sys
import argparse
from numpy import save, load

sys.path.insert(0, '../src')
from hamiltonian import *
from misc import *
from haldane_pseudopotential import *


parser = argparse.ArgumentParser(description='Exact diagonalize of spherical fractional quantum hall system.')

parser.add_argument('--Ne', type=int, required=True,
                    help='(int) Number of electrons')
parser.add_argument('--m', type=int, required=True,
                    help='(int) Filling factor v=1/m fulfilling conditions (0<v<1). If v exceeds 1 (such as the 5/2 pfaffian state), enter the highest partial-filled LL filling factor.')
parser.add_argument('--wf', type=str, required=True,
                    help='(float) L: Laughlin state (S=unity), MR-Pf: Moore-Read Pfaffian state (S=3), aPf: anti-Pfaffian state (S=1), PH-Pf: particle-hole conjugation (S=-1)')
parser.add_argument('--nLL',  nargs='?', const=0, type=int,
                    help='(int, opt) nth Landau level. Default value nLL=0 if not specified')
parser.add_argument('--L', nargs='?', const=1, type=int,
                    help='(int, opt) Angular momentum. Default value L=0 if not specified.')
parser.add_argument('--ppm', type=str, required=True,
                    help='(str, opt) Filename for input two-body interaction matrix.')
parser.add_argument('--type',  type=str, required=True,
                    help='(str, opt) Name for interaction. Has no effect on calculations. Only provides labeling for output file name.')
parser.add_argument('--getE', nargs='?', default=False, type=bool,
                    help='(bool, opt) Save eigenenergy to file. Default set to Fasle.')
parser.add_argument('--getL2', nargs='?', default=False, type=bool,
                    help='(bool, opt) Generate L2 matrix and save to file. Default set to False.')
parser.add_argument('--getHamil', nargs='?', default=False, type=bool,
                    help='(bool, opt) Save full Hamiltonian to file. Default set to Fasle.')
args = parser.parse_args()

wf_states={'L':args.m, 'MR-Pf':3, 'aPf':1, 'PH-Pf':-1}
for state in wf_states:
    if args.wf==state: 
        S=wf_states[state]

N=args.Ne
nLL=args.nLL
Nphi=args.m*N-S

system=system(N,Nphi,0,args.getL2)

pp_matrix=load(args.ppm)

hamiltonian = spherical_system(system, pp_matrix)
energies, eigenstates=get_eig(hamiltonian)

if args.getL2==True:
    save(f'{args.wf}{args.Ne}_{args.type}_L2', system.L2.A)

if args.getE==True:
    save(f'{args.wf}{args.Ne}_{args.type}_energies', energies)

if args.getHamil==True:
    save(f'{args.wf}{args.Ne}_{args.type}_Hamil', hamiltonian.h.A)

save(f'{args.wf}{args.Ne}_{args.type}_eigenstates', eigenstates)



