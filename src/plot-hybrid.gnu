set terminal pngcairo size 800,600 enhanced font "Arial,12"
set output "runtime-hybrid.png"

set title "Runtime vs. Processes and Threads" font ",14"
set xlabel "Number of Processes" font ",12"
set ylabel "Number of Threads per Process" font ",12"
set zlabel "Runtime (s)" font ",12"
set grid

set dgrid3d 20,20,1  # Smooth the surface
set hidden3d         # Enable hidden surface removal
set ticslevel 0      # Set the z-axis to align with the base

splot "output.txt" using 1:2:3 with linespoints palette title "Runtime"
