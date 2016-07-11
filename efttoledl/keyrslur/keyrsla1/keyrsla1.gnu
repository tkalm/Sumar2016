reset

# Output
set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'plotk1.eps'
set size 2,1

# Axes
set xrange [0:100]
set yrange [0:1.02]
set xlabel "Time [ps]"
set ylabel "Occupation"

# Header
set title "Time evolution of the occupation of states \n \
{Initalial state |0>, {/Symbol k} = 0.2E_0, \n \
Energies: E_0 = {/=12 @^{/Symbol=10 -}{/=14 h}}{/Symbol w}_0 = 1.0 meV \
and  {/=12 @^{/Symbol=10 -}{/=14 h}}{/Symbol W} = 0.4 meV"

# Key
set key right center 
set key font ",18"
set key spacing "1.8"

# Line colors
set style line 2 lc rgb "#D4A006" lt 1 lw 1.5
set style line 3 lc rgb "#D4A006" lt 1 lw 1.5
set style line 4 lc rgb "#D4A006" lt 1 lw 1.5
set style line 5 lc rgb "#D4A006" lt 1 lw 1.5

# Plot 
plot for [col=2:5] "nit.dtx" using 1:col with lines ls col title columnheader,\
