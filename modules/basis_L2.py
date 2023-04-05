#!/usr/bin/env python

import sys
import argparse
from numpy import save, load

sys.path.insert(0, '../src')
from misc import *

parser = argparse.ArgumentParser(description='Generates energetics of specified FQHE system in terms of angular momentum. Good for visual representation of states for spherical systems.')


parser.add_argument('--label', type=str, required=True,
                    help='(str) Label for filename. No impact on system.')
parser.add_argument('--H', type=str, required=True,
                    help='(str) Filename for Hamiltonian of system.')
parser.add_argument('--L2', type=str, required=True,
                    help='(str) Filename for L2 matrix.')
args = parser.parse_args()

hamil=load(args.H)
L2=load(args.L2)

L,ene = L_spectrum(L2,hamil)

print('==========================================================')
print(f'Energy eigenstates organized by angular momentum L saved to file {args.label}_LBasis.npy')
print('')
save(f'{args.label}_LBasis', (L,ene))