#!/usr/bin/env gnuplot --persist

set terminal qt 0 size 500,400 font ",20"

set key left

set mxtics
set mytics

set xlabel "Temperature [K]"
set ylabel "Pressure [bar]"

plot "output/env-2ph-PT_1.dat" u 4:5 w l lc "blue" ,\
     "output/env-2ph-PT_2.dat" u 4:5 w l lc "black" ,\

pause mouse close