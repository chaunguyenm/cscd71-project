#!/bin/bash
#SBATCH --nodes=4
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --job-name stft_mpi_job
#SBATCH --output=stft_mpi_output_%j.txt
#SBATCH --mail-type=FAIL

cd $SLURM_SUBMIT_DIR

# load required modules
source teachsetup

# run scaling analysis
for i in {1..4}; do
    { time -p mpirun -n "$i" "$@"; } 2> "output_$i.txt"
done

# prepare output
outputs=$(find . -maxdepth 1 -name "output_*.txt" -print | sort -t_ -k2 -n)

for output in $outputs; do
    nthreads=$(echo $output | sed 's/.*_\([0-9]\+\)\.txt/\1/')
    runtime=$(grep "real" "$output" | sed 's/real //')

    if [ "$nthreads" -eq 1 ]; then
        runtime1=$runtime
        echo "$nthreads $runtime 1" > "output.txt"
    else
        speedup=$(echo "$runtime1 / $runtime" | bc -l)
        echo "$nthreads $runtime $speedup" >> "output.txt"
    fi
done

# plot
module load gnuplot
gnuplot -e "dir='.'" plot.gnu

