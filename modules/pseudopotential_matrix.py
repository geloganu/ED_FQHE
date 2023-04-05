#!/usr/bin/env python
import sys
import argparse
from numpy import save

sys.path.insert(0, '../src')
from haldane_pseudopotential import *

parser = argparse.ArgumentParser(description='Reads pseudopotential Vm and generates the two-body interaction matrix in spherical geometry.')

parser.add_argument('--Nphi', type=int, required=True,
                    help='(int) Number of flux quanta on the surface of the sphere.')
parser.add_argument('--nLL',  nargs='?', const=0, type=int,
                    help='(int, opt) nth Landau level. Default value nLL=0 if not specified')
parser.add_argument('--interaction', type=str, required=True,
                    help='(str) C:Coulomb, T:Trial')
parser.add_argument('--input',  nargs='?', const=None, 
                    help='(str, opt) If interaction type T, require input pseudopotential file')
args = parser.parse_args()

if args.interaction=='C':
    pp = haldane_pseudopotential(args.Nphi/2, LLn = args.nLL)
elif args.interaction=='T':
    Vm=np.loadtxt(args.input)
    pp = haldane_pseudopotential(args.Nphi/2, LLn = args.nLL,custom=Vm)

save(f'Nphi{args.Nphi}_{args.interaction}', pp.pp_matrix)