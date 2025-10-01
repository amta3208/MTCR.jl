#!/bin/bash

quiet=0

while getopts 'qn:' flag; do
    case "${flag}" in
        n) 
            nprocs=${OPTARG}
            ;;
        q) 
            quiet=1
            ;;
    esac
done

if [[ -d ./input ]]; then
    if [[ -f ./input/prob_setup.inp ]]; then
        echo "Starting MTCR run: $(date)" 
        if [[ -d ./output ]]; then
            rm -r ./output
        fi
        mkdir ./output
        mkdir ./output/states
        mkdir ./output/sources
    else
        echo "Could not find ./input/prob_setup.inp"; exit 1
    fi
else
    echo "Could not find ./input"; exit 1
fi 

if [ $quiet -eq 1 ]; then
    mpirun -n $nprocs "$MTCR_INSTALL_DIR/mtcr.x" > printout.out
else
    echo ""
    mpirun -n $nprocs "$MTCR_INSTALL_DIR/mtcr.x" | tee printout.out
    echo ""
fi

echo "Finished MTCR run: $(date)"
