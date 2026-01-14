#!/bin/bash

# Clean up temporary files and directories created during the assignment
echo "Cleaning up make artifacts..."
make clean
echo "All .dat files have been removed."
rm -f experiments/*.dat
#rm -f run_experiments_*.out
#echo "All run_experiments output files have been removed."