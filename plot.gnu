#!/usr/bin/env gnuplot --persist

set terminal qt 0 size 500,400 font ",20"

set key left

set mxtics
set mytics

set xlabel "Temperature [K]"
set ylabel "Pressure [bar]"

plot "ENV2_OUT_1" index 0 using 1:2 with lines linecolor rgb "black" title "inc: Vapor", \
     "ENV2_OUT_1" index 1 using 1:2 with lines linecolor rgb "blue" title "inc: Liquid", \
     "ENV2_OUT_1" index 2 using 1:2 with points linestyle 7 lc rgb "black" title "CP"

pause mouse close
