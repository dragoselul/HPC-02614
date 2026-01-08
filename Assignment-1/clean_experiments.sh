#!/bin/bash

# Script to clean up generated .dat files from matrix multiplication experiments
# Usage: ./clean_experiments.sh
# This script removes all .dat files created by the run_experiments.sh script
rm -f *.dat
echo "All .dat files have been removed."
rm -f run_experiments_*.out
echo "All run_experiments output files have been removed."