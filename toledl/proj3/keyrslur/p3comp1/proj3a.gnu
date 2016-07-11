reset

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'graphproj3.eps'

set xrange [0:100]
set yrange [0:1.03]

set xlabel "Time [ps]"
set ylabel "Occupation"
set title "Time evolution of occupation: comparison of direct calculation and time integration \
in a system with damping strength {/Symbol k} = 0.1E_0\
\n Initalial state |0>, Energies: E_0 = {/=12 @^{/Symbol=10 -}{/=14\
h}}{/Symbol w}_0 = 1.0 meV and  {/=12 @^{/Symbol=10 -}{/=14 h}}{/Symbol W} = 0.4 meV"
set size 2,1
set key right center 
set key font ",18"
set key spacing "1.8"

plot for [col=2:4] "nit.dtx" using 1:col with lines lt rgb "#E8D900" title "Direct",\
     "aState0.dtx" with lines lt rgb "#006400" title "Time integrated",\
     "aState1.dtx" with lines lt rgb "#006400" title "Time integrated",\
     "aState2.dtx" with lines lt rgb "#006400" title "Time integrated",\
