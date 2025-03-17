#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --job-name stft_hybrid_job
#SBATCH --output=stft_hybrid_output_%j.txt
#SBATCH --mail-type=FAIL

cd $SLURM_SUBMIT_DIR

# load required modules
source teachsetup

# run scaling analysis
for nprocs in {1..8}; do
    for nthreads in {1..8}; do
        export OMP_NUM_THREADS=$nthreads
        { time -p mpirun -n "$nprocs" "$@"; } 2> "output_${nprocs}_${nthreads}.txt"
    done
done

rm -f output.txt

# prepare output
outputs=$(find . -maxdepth 1 -name "output_*_*.txt" -print | sort -t_ -k2,2n -k3,3n)

for output in $outputs; do
    nthreads=$(echo $output | sed 's/.*_\([0-9]\+\)\.txt/\1/')
    nprocs=$(echo $output | sed 's|.*/output_\([0-9]\+\)_.*|\1|') 
    runtime=$(grep "real" "$output" | sed 's/real //')
    echo "$nprocs $nthreads $runtime" >> "output.txt"
done

# plot
module load gnuplot
gnuplot plot-hybrid.gnu

