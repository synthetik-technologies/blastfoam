#! /bin/bash
gnuplot -p << EOF
set grid
set title 'Displacement of the Flap Tip'
set xlabel 'Time [s]'
set ylabel 'X-Displacement [m]'
# set style line 1 lt 2 lw 6
# set style line 1 lt 2 lw 6
# set linestyle  2 lt 2 lc 1 # red-dashed
set linestyle  1 lt 2 lc 1 # red-dashed
plot "precice-Fluid-watchpoint-point1.log" using 1:8 title 'Top displacemement' with lines
EOF

