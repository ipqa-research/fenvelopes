#!/usr/bin/gnuplot -persist

set key left
set mytics 4
set mxtics 5

set xlabel "x"
set ylabel "P [bar]"

set xrange [0:1]
set xrange [0:1]
set yrange [0:2000]

dew_2ph = 1
bub_2ph = 2
dew_3ph = 3
bub_3ph = 4

set style line dew_2ph lc rgb "blue"
set style line bub_2ph lc rgb "black"

set style line dew_3ph lc rgb "blue" dashtype 2
set style line bub_3ph lc rgb "black"  dashtype 2


plot "./output/env-2ph-PX_1.dat" u 4:5 w l ls bub_2ph t "Bubble Line", \
     "./output/env-2ph-PX_2.dat" u 4:5 w l ls dew_2ph t "Dew Line", \
     "./output/env-3ph-PX_3.dat" u 4:5 w l ls bub_3ph t "3ph-vapor", \
     "./output/env-3ph-PX_4.dat" u 4:5 w l ls dew_3ph t "3ph-liquid", \
     "./output/env-3ph-PX_3.dat" index "critical" u 1:2 w p t "" pt 7 lc rgb "black", \
     "./output/env-3ph-PX_4.dat" index "critical" u 1:2 w p t "" pt 7 lc rgb "black"
pause mouse close
