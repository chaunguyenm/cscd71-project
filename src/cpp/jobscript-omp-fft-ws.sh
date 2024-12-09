#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --job-name stft_omp_job
#SBATCH --output=stft_omp_output_%j.txt
#SBATCH --mail-type=FAIL

cd $SLURM_SUBMIT_DIR

# load required modules
source teachsetup

num_samples=131072

num_threads=1
export OMP_NUM_THREADS=$num_threads
window_size=256
args="$num_samples $window_size"
{ time -p "$@" $args; } 2> "output_$num_threads.txt"

num_threads=2
export OMP_NUM_THREADS=$num_threads
window_size=512
args="$num_samples $window_size"
{ time -p "$@" $args; } 2> "output_$num_threads.txt"

num_threads=8
export OMP_NUM_THREADS=$num_threads
window_size=1024
args="$num_samples $window_size"
{ time -p "$@" $args; } 2> "output_$num_threads.txt"

num_threads=24
export OMP_NUM_THREADS=$num_threads
window_size=2048
args="$num_samples $window_size"
{ time -p "$@" $args; } 2> "output_$num_threads.txt"

num_threads=64
export OMP_NUM_THREADS=$num_threads
window_size=4096
args="$num_samples $window_size"
{ time -p "$@" $args; } 2> "output_$num_threads.txt"

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
