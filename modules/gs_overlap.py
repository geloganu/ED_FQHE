#!/usr/bin/env python

import sys
import argparse
from numpy import load, transpose

sys.path.insert(0, '../src')
from misc import *

parser = argparse.ArgumentParser(description='Overlap of groundstate eigenstates. Vector overlap (dot product).')

parser.add_argument('--vec1', type=str, required=True,
                    help='(str) Filename for first vector')
parser.add_argument('--vec2', type=str, required=True,
                    help='(str) Filename for second vector')
args = parser.parse_args()

vec1=load(args.vec1)
vec2=load(args.vec2)

overlap=overlap_Vectors(np.transpose(vec1)[0],np.transpose(vec2)[0])