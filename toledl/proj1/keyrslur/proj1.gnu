reset
set xrange [0:1]
set yrange [0:6]
set xlabel "{/Symbol l}"
set ylabel "E_n/h{/Symbol n}"
set title "Energy spectrum of H = H_0 + {/Symbol l}V"
set size 1.8,1 

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'graphproj1.eps'

plot "State0.dtx" with lines lt rgb "red" notitle,\
     "State1.dtx" with lines lt rgb "#006400" notitle,\
     "State2.dtx" with lines lt rgb "blue" notitle,\
     "State3.dtx" with lines lt rgb "#FF1493" notitle,\
     "State4.dtx" with lines lt rgb "cyan" notitle,\
     "State5.dtx" with lines lt rgb "#D2691E" notitle,\
