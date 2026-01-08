#!/bin/bash
# Script to generate .dat files for different matrix multiplication algorithms
# Similar to run_tunelab.sh from day-2/lab
#BSUB -J run_experiments
#BSUB -o run_experiments_%J.out
#BSUB -q hpcintro
#BSUB -n 1
#BSUB -R "rusage[mem=2048]"
#BSUB -W 15
#BSUB -R "span[hosts=1] affinity[socket(1)]"

EXECUTABLE=matmult_c.gcc
SIZES="50 100 150 200 300 400 500 700 1000 1500 2000"
PERMS="nat lib mnk mkn kmn knm nmk nkm"
BLKSIZES="8 16 32 64 128 256"

# Clean old data files
rm -f *.dat

# Run permutation experiments
for PERM in $PERMS; do
    echo "Running $PERM..."
    for S in $SIZES; do
        ./$EXECUTABLE $PERM $S $S $S >> ${PERM}.dat
    done
done

# Run block size experiments
for BS in $BLKSIZES; do
    echo "Running blk with bs=$BS..."
    for S in $SIZES; do
        ./$EXECUTABLE blk $S $S $S $BS >> blk_${BS}.dat
    done
done

echo "Done! Data files generated."

