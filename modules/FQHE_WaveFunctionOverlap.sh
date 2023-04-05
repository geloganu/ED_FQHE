#!/usr/bin/env bash

if [ "$1" == "-h" ]; then
  echo -e "\nThis programs calculates the overlap between the exact Coulomb and trial wavefunction eigenstates on the Haldane sphere."
  echo "==================================================================="
  echo "IMPORTANT:"
  echo "The trial pseudopotentials to generate the trial Hamiltonian should be named using the following convention:Nphi[Nphi]_pp.txt where [Nphi] is the total number of flux quanta on the sphere."
  echo "==================================================================="
  echo "ARGS:"
  echo "    -e   Ne:   (int) Number of electrons"
  echo "    -p   Nphi: (int) Number of flux quanta on the surface of the Haldane sphere"
  echo "    -n   nLL:  (int) nth Landau level"
  echo "    -m   m:    (int) Filling factor v=1/m"
  echo "    -w   wf:   (str) L: Laughlin, MR-Pf: Moore-Read Pf, aPf: anti-Pf, PH-Pf: PH conjugation"
  echo "==================================================================="
  echo "EXAMPLE:"
  echo "    *Create trial pseudopotential file names 'Nphi21_pp.txt' wit the following contents: 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0" 
  echo "    *To diagonalize the Laughlin state for Ne=8 electrons: PATH/TO/RunFQHE.sh -e 8 -p 21 -n 0 -m 3 -w L"
  echo " "
  exit 0
fi

while getopts e:p:n:m:w: flag;
do
  case "${flag}" in
    e)
      Ne=${OPTARG}
      ;;
    p)
      Nphi=${OPTARG}
      ;;
    n)
      nLL=${OPTARG}
      ;;
    m)
      m=${OPTARG}
      ;;
    w)
      wf=${OPTARG}
      ;;
  esac
done

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

#exact diagonalize exact wavefunction
python $SCRIPT_DIR/pseudopotential_matrix.py --Nphi ${Nphi} --nLL ${nLL} --interaction C
python $SCRIPT_DIR/fqhe_ed.py --Ne ${Ne} --m ${m} --wf ${wf} --ppm Nphi${Nphi}_C.npy --type Coulomb 

#exact diagonalize trial wavefunction
python $SCRIPT_DIR/pseudopotential_matrix.py --Nphi ${Nphi} --nLL ${nLL} --interaction T --input Nphi${Nphi}_pp.txt
python $SCRIPT_DIR/fqhe_ed.py --Ne ${Ne} --m ${m} --wf ${wf} --ppm Nphi${Nphi}_T.npy --type Trial

#yield overlap
python $SCRIPT_DIR/gs_overlap.py --vec1 ${wf}${Ne}_Coulomb_eigenstates.npy --vec2 ${wf}${Ne}_Trial_eigenstates.npy
