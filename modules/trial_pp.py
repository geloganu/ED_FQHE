import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Generate trial pseudopotential Vm values for specified m.')

parser.add_argument('--Nphi', type=int, required=True,
                    help='(int) Number of flux quanta on the surface of the sphere.')
parser.add_argument('--Vm', nargs='*', type=str, required=True,
                    help='(int) Input values for Vm.')
args = parser.parse_args()

input= [num for num in args.Vm[0].split(":")]

if len(input)%2 != 0:
    raise ValueError('Missing input in Vm')

pos=np.asarray(input[::2],dtype=int) #list denoting position of custom Vm
Vm=input[1::2] #list contaiining custom Vm

trial_pp=np.zeros(args.Nphi+1)
trial_pp[pos]=Vm

print('Trial Pseudopotentials:', trial_pp, f'saved in file Nphi{args.Nphi}_pp.txt')

np.savetxt(f'Nphi{args.Nphi}_pp.txt',trial_pp.reshape(1, trial_pp.shape[0]))