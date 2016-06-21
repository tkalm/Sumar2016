reset
#set xrange [0:67]
#set yrange [0:1.1]
set output 'graphproj2b1.eps'
set xlabel "Time [ps]"
set ylabel "<x>"
set title "Time evolution of the expectation value of x"
set size 2,1

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14

plot "expx.dtx" with lines lt rgb "green" notitle
