#!/usr/bin/env bash

if [ "$1" == "-h" ]; then
  echo -e "\nThis program uses exact diagonalization to simulate and retrieve the energetics of the specified FQHE system."
  echo "==================================================================="
  echo "ARGS:"
  echo "    -e   Ne:    (int) Number of electrons"
  echo "    -p   Nphi:  (int) Number of flux quanta on the surface of the Haldane sphere"
  echo "    -n   nLL:   (int) nth Landau level"
  echo "    -m   m:     (int) Filling factor v=1/m"
  echo "    -w   wf:    (str) L: Laughlin, MR-Pf: Moore-Read Pf, aPf: anti-Pf, PH-Pf: PH conjugation"
  echo "==================================================================="
  echo "EXAMPLE:"
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
SECONDS=0

#exact diagonalize exact wavefunction
python $SCRIPT_DIR/pseudopotential_matrix.py --Nphi ${Nphi} --nLL ${nLL} --interaction C
python $SCRIPT_DIR/fqhe_ed.py --Ne ${Ne} --m ${m} --wf ${wf} --ppm Nphi${Nphi}_C.npy --type Coulomb --getL2 True --getHamil True
python $SCRIPT_DIR/basis_L2.py --label ${wf}${Ne} --H ${wf}${Ne}_Coulomb_Hamil.npy --L2 ${wf}${Ne}_Coulomb_L2.npy

echo $SECONDS