#!/usr/bin/gnuplot -persist

set key left
set mytics 4
set mxtics 5

set xlabel "x"
set ylabel "P [bar]"

plot "./output/env-2ph-PX_1.dat" u 4:5 w lp t "Bubble Line", \
     "./output/env-2ph-PX_2.dat" u 4:5 w lp t "Dew Line", \
     "./output/env-3ph-PX_3.dat" u 4:5 w lp t "3ph-vapor", \
     "./output/env-3ph-PX_4.dat" u 4:5 w lp t "3ph-liquid", \
     "./output/env-3ph-PX_3.dat" index "critical" u 1:2 w p t "" pt 7 lc rgb "black", \
     "./output/env-3ph-PX_4.dat" index "critical" u 1:2 w p t "" pt 7 lc rgb "black"

pause mouse close
