#!/bin/bash

set -e
source "$MTCR_INSTALL_DIR/tools/postprocess.sh"

# Initializations
calling_folder=$(pwd)
mtcr_driver_script="$MTCR_INSTALL_DIR/tools/run_mtcr.sh"

# Default variables
nprocs=1
overwrite=$true # Overwrite .mat file if one is found

# Read command line arguments if found
while getopts n: flag
do
    case "${flag}" in
        n) nprocs=${OPTARG};;
    esac
done
for arg in "$@"; do
    if [ "$arg" = "-no_overwrite" ]; then 
        overwrite=$false
    fi
done

# Loop through directories and run each case
for folder in */; do
    case_folder=${folder%/}
    # run the case  
    echo "Running case: $case_folder"
    cd "$case_folder"
    bash "$mtcr_driver_script" -n $nprocs
    # run postprocessor if: no mat file is found -or- if overwrite is true
    mat_file_name="$case_folder.mat" 
    if  ! [ -f "$mat_file_name" ] || [ ! $overwrite ]; then
        rm -f "$mat_file_name"
        postprocess_mtcr
    fi
    cd "$calling_folder"
done