reset
set xlabel "Time [ps]"
set ylabel "Occupation"
set title "Occupation of the lowest states of H"
set size square

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'OccStatPlot.eps'

plot "state0" with lines notitle,\
     "state1" with lines notitle,\
     "state2" with lines notitle,\
     "state3" with lines notitle,\
     "state4" with lines notitle,\
