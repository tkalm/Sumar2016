reset
set xlabel "Time [ps]"
set ylabel "Occupation"
set title "The Occupation of the five lowest states"
set size square

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'proj2plot.eps'

plot "State0.dtx" with lines title "|0>",\
     "State1.dtx" with lines title "|1>",\
     "State2.dtx" with lines title "|2>",\
     "State3.dtx" with lines title "|3>",\
     "State4.dtx" with lines title "|4>",\
