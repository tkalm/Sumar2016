reset
set xrange [0:100]
set yrange [0:1.1]
set xlabel "Time [ps]"
set ylabel "Occupation"
set title "Time evolution of the occupation of the five lowest states of H = H_{0} + H'(t) - Dissipation"
set size 2,1

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'graphproj2c.eps'

plot "trace.dtx" with lines lt rgb "black" title "Trace",\
     "State0.dtx" with lines lt rgb "red" title "|0>",\
     "State1.dtx" with lines lt rgb "green" title "|1>",\
     "State2.dtx" with lines lt rgb "blue" title "|2>",\
     "State3.dtx" with lines lt rgb "purple" title "|3>",\
     "State4.dtx" with lines lt rgb "cyan" title "|4>"
