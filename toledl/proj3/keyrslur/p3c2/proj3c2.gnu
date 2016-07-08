reset

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'graphproj3c2.eps'

set xrange [0:105]
set yrange [0:1.03]
set xlabel "Time [ps]"
set ylabel "Occupation"
set title "Time evolution of the occupation of states {/Symbol m}_n in a system with damping strength {/Symbol k}\
\n Initalial state |1>, Energies: E_0 = {/=12 @^{/Symbol=10 -}{/=14\
h}}{/Symbol w}_0 = 1.0 meV and  {/=12 @^{/Symbol=10 -}{/=14 h}}{/Symbol W} = 0.4 meV"
set size 2,1
set key right center 
set key font ",18"
set key spacing "1.8"

plot "aState0.dtx" with lines lt rgb "#006400" title "{/Symbol k} = 0.1E_0    {/Symbol m}_0",\
     "aState1.dtx" with lines lt rgb "#006400" title "{/Symbol k} = 0.1E_0    {/Symbol m}_1",\
     "aState2.dtx" with lines lt rgb "#006400" title "{/Symbol k} = 0.1E_0    {/Symbol m}_2",\
     "aState3.dtx" with lines lt rgb "#006400" title "{/Symbol k} = 0.1E_0    {/Symbol m}_3",\
     "bState0.dtx" with lines lt rgb "#FF8C00" title "{/Symbol k} = 0.05E_0  {/Symbol m}_0",\
     "bState1.dtx" with lines lt rgb "#FF8C00" title "{/Symbol k} = 0.05E_0  {/Symbol m}_1",\
     "bState2.dtx" with lines lt rgb "#FF8C00" title "{/Symbol k} = 0.05E_0  {/Symbol m}_2",\
     "bState3.dtx" with lines lt rgb "#FF8C00" title "{/Symbol k} = 0.05E_0  {/Symbol m}_3",\
     "steadystate.dtx" using ($1*0+100):2 with points pt 13 ps 1.8 title "Steady states",\
