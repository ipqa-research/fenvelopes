#!/usr/bin/gnuplot -persist

set terminal qt font ",16"
set key font ",13"

set key left
set mytics 4
set mxtics 5

set xlabel "{/Symbol a}"
set ylabel "P [bar]"

set xrange [0:1]
set xrange [0:1]
set yrange [0:2000]

dew_2ph = 1
bub_2ph = 2
liq_2ph = 3

dew_3ph = 4
bub_3ph = 5
liq_3ph = 6

set style line dew_2ph lc rgb "blue"  lw 1.5
set style line bub_2ph lc rgb "black" lw 1.5
set style line liq_2ph lc rgb "red"   lw 1.5

set style line dew_3ph lc rgb "blue"  dt 7 lw 1.5
set style line bub_3ph lc rgb "black" dt 7 lw 1.5
set style line liq_3ph lc rgb "red"   dt 7 lw 1.5

pxs_2 = system("ls -d ./fenvelopes_output/env-2ph-PX* | xargs")
pxs_3 = system("ls -d ./fenvelopes_output/env-2ph-PX* | xargs")

plot "fenvelopes_output/env-2ph-PX_1.dat" u 4:5 w lp ls bub_2ph t "2ph-bub",\
     "fenvelopes_output/env-2ph-PX_2.dat" u 4:5 w lp ls dew_2ph t "2ph-dew", \
     "fenvelopes_output/env-2ph-PX_3.dat" u 4:5 w lp ls dew_2ph t "2ph-dew", \
     "fenvelopes_output/env-2ph-PX_4.dat" u 4:5 w lp ls dew_2ph t "2ph-dew", \
     "fenvelopes_output/env-3ph-PX_3.dat" u 4:5 w lp ls bub_3ph t "3ph-bub", \
     "fenvelopes_output/env-3ph-PX_4.dat" u 4:5 w lp ls dew_3ph t "3ph-dew", \
     "fenvelopes_output/env-3ph-PX_5.dat" u 4:5 w lp ls bub_3ph t "3ph-bub", \
     "fenvelopes_output/env-3ph-PX_6.dat" u 4:5 w lp ls dew_3ph t "3ph-dew", \

pause mouse close
