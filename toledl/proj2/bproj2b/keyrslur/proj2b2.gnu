reset
#set xrange [0:67]
#set yrange [0:1.1]
set output 'graphproj2b2.eps'
set xlabel "Time [ps]"
set ylabel "<{x^2}>"
set title "Time evolution of the expectation value of {x^2}"
set size 2,1

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14

plot "expxx.dtx" with lines lt rgb "red" notitle
