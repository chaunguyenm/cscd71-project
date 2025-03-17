if (!exists("dir")) dir='.'

set terminal pngcairo
set output dir."/runtime-speedup.png"
set multiplot layout 2, 1 title "Run time and Speedup by core number"

set title "Run time by core number"
set xlabel "Core count"
set ylabel "Run time (s)"
plot dir."/output.txt" using 1:2 title "Run time"  with lines

set title "Speedup by core number"
set xlabel "Core count"
set ylabel "Speedup"
plot dir."/output.txt" using 1:3 title "Speedup" with lines
