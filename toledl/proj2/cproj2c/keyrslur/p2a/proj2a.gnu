reset

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'graphproj2a.eps'

set xrange [0:100]
set yrange [0:1.2]
set xlabel "Time [ps]"
set ylabel "Occupation"
set title "Time evolution of the occupation of states {/Symbol m}_n in a system with damping strength {/Symbol k} = 0\
 eV \n Initalial state |0>, Energies: E_0 = {/=12 @^{/Symbol=10 -}{/=14\
h}}{/Symbol w}_0 = {/=12 @^{/Symbol=10 -}{/=14 h}}{/Symbol W} = 1.0 meV"

set size 2,1

plot "trace.dtx" with lines lt rgb "black" title "Trace",\
     "State0.dtx" with lines lt rgb "red" title "{/Symbol m}_0",\
     "State1.dtx" with lines lt rgb "green" title "{/Symbol m}_1",\
     "State2.dtx" with lines lt rgb "blue" title "{/Symbol m}_2",\
     "State3.dtx" with lines lt rgb "purple" title "{/Symbol m}_3",\
     "State4.dtx" with lines lt rgb "cyan" title "{/Symbol m}_4"
