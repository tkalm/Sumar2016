reset

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'graphproj2c3.eps'

set xrange [0:100]
set yrange [0:1.03]
set xlabel "Time [ps]"
set ylabel "Occupation"
set title "Time evolution of the occupation of states {/Symbol m}_n in a system with damping strength {/Symbol k}\
 = 0.1E_0 \n Initalial state |2>, Energies: E_0 = {/=12 @^{/Symbol=10 -}{/=14\
h}}{/Symbol w}_0 = 1.0 meV and  {/=12 @^{/Symbol=10 -}{/=14 h}}{/Symbol W} = 0.6 meV"
set size 2,1
set key right center 
set key font ",18"
set key spacing "1.8"

plot "trace.dtx" with lines lt rgb "#2C5A8F" title "Trace",\
     "State0.dtx" with lines lt rgb "red" title "{/Symbol m}_0",\
     "State1.dtx" with lines lt rgb "green" title "{/Symbol m}_1",\
     "State2.dtx" with lines lt rgb "blue" title "{/Symbol m}_2",\
     "State3.dtx" with lines lt rgb "purple" title "{/Symbol m}_3",\
     "State4.dtx" with lines lt rgb "cyan" title "{/Symbol m}_4",\
