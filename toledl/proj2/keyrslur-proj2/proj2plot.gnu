reset
set xlabel "Time [ps]"
set ylabel "Occupation"
set title "The Occupation of the five lowest states"
set size square

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'proj2plot.eps'

plot "state0" with lines title "|0>",\
     "state1" with lines title "|1>",\
     "state2" with lines title "|2>",\
     "state3" with lines title "|3>",\
     "state4" with lines title "|4>",\
