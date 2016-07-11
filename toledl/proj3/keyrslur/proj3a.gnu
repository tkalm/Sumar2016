reset

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'graphproj3.eps'

set xrange [0:100]
set yrange [0:1.03]

#set logscale y
set xlabel "Time [ps]"
set ylabel "Occupation"
set title "Time evolution"
set size 2,1
set key right center 
set key font ",18"
set key spacing "1.8"

plot for [col=2:16] "nit.dtx" using 1:col with lines title columnheader
