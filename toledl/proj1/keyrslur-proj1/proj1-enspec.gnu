reset
set xrange [0:1]
set yrange [0:6]
set xlabel "{/Symbol l}"
set ylabel "E_n/h{/Symbol n}"
set title "Energy spectrum of H = H_0 + {/Symbol l}V"
set size square

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'proj1-enspec.eps'

plot "state1" with lines notitle,\
     "state2" with lines notitle,\
     "state3" with lines notitle,\
     "state4" with lines notitle,\
     "state5" with lines notitle,\
     "state6" with lines notitle\
