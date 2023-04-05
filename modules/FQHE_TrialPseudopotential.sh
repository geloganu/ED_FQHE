#!/usr/bin/env bash

if [ "$1" == "-h" ]; then
  echo -e "\nThis programs generates a file that contains the trial pseudopotentials."
  echo "==================================================================="
  echo "ARGS:"
  echo "    -p   Nphi: (int) Number of flux quanta on the surface of the Haldane sphere"
  echo "    -v   Vm:   (int) Takes multiples of 2 inputs. The first value is the position of the custom"
  echo "                     Vm value followed by the value of Vm. Repeat for other custom values"
  echo "==================================================================="
  echo "EXAMPLE:"
  echo "    *For the Laughlin state with Nphi=15, use: ../modules/FQHE_TrialPseudopotential.sh -p 15 -v 0:1:1:1" 
  echo "This generates 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0."
  echo ""
  exit 0
fi

while getopts p:v: flag;
do
  case "${flag}" in
    p)
      Nphi=${OPTARG}
      ;;
    v)
      Vm=${OPTARG}
      ;;
    
  esac
done

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
SECONDS

#exact diagonalize exact wavefunction
python $SCRIPT_DIR/trial_pp.py --Nphi ${Nphi} --Vm ${Vm}

echo $SECONDS