reset

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'graphproj2d1.eps'

#set xrange [0:100]
#set yrange [0:1.03]
set xlabel "Time [ps]"
set ylabel "x"
set title "Expectation value of x and x^2 in a system with damping strength {/Symbol k}\
 = 0 eV \n Initalial state |0>, Energies: E_0 = {/=12 @^{/Symbol=10 -}{/=14\
h}}{/Symbol w}_0 = {/=12 @^{/Symbol=10 -}{/=14 h}}{/Symbol W} = 1.0 meV"
set size 2,1
#set key right center 
#set key font ",18"
#set key spacing "1.8"

plot "expx.dtx" with lines lt rgb "#006400" title "Expectation of x",\
     "expxx.dtx" with lines lt rgb "#FF8C00" title "Expectation of x^2",\
