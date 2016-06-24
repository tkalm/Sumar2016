reset

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'graphproj2c.eps'

set xrange [0:100]
set yrange [0:1.2]
set xlabel "Time [ps]"
set ylabel "Occupation"
set title "Time evolution of the occupation of states with damping strength  {/Symbol k}"
set size 2,1

plot "trace.dtx" with lines lt rgb "black" title "Trace",\
     "State0.dtx" with lines lt rgb "red" title "E_0 = h{/Symbol n} = 1.0 meV, {/Symbol k} = 0.1E_0, , R=1.0",\
     "State1.dtx" with lines lt rgb "green" title "|1>",\
     "State2.dtx" with lines lt rgb "blue" title "|2>",\
     "State3.dtx" with lines lt rgb "purple" title "|3>",\
     "State4.dtx" with lines lt rgb "cyan" title "|4>"
