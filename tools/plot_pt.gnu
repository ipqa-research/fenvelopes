#!/usr/bin/env gnuplot --persist

set terminal qt 0 size 500,400 font ",20"

set key left

set mxtics
set mytics

set xlabel "Temperature [K]"
set ylabel "Pressure [bar]"

plot "fenvelopes_output/env-2ph-PT_1.dat" u 4:5 w lp,\
     "fenvelopes_output/env-2ph-PT_2.dat" u 4:5 w lp,\
     "fenvelopes_output/env-2ph-PT_3.dat" u 4:5 w lp,\
     "fenvelopes_output/env-3ph-PT_4.dat" u 4:5 w lp,\
     "fenvelopes_output/env-3ph-PT_5.dat" u 4:5 w lp,\
     "fenvelopes_output/env-3ph-PT_6.dat" u 4:5 w lp,\
     "fenvelopes_output/env-3ph-PT_7.dat" u 4:5 w lp,\

#plot "fenvelopes_output/env-2ph-PT_1.dat" u 4:5 w l lc "black" ,\
#     "fenvelopes_output/env-2ph-PT_2.dat" u 4:5 w l lc "blue" ,\
#     "fenvelopes_output/env-2ph-PT_3.dat" u 4:5 w l lc "blue" ,\
#     "fenvelopes_output/env-2ph-PT_4.dat" u 4:5 w l lc "blue" ,\
#     "fenvelopes_output/env-3ph-PT_3.dat" u 4:5 w l lc "black" ,\
#     "fenvelopes_output/env-3ph-PT_4.dat" u 4:5 w l lc "blue" ,\
#     "fenvelopes_output/env-3ph-PT_5.dat" u 4:5 w l lc "black" ,\
#     "fenvelopes_output/env-3ph-PT_6.dat" u 4:5 w l lc "blue" ,\

pause mouse close
