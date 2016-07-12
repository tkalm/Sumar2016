reset

# Output
set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'plotnxk2.eps'
set size 2,1

# Axes
set xrange [0:100]
set yrange [0:1.02]
set xlabel "Time [ps]"
set ylabel "Occupation"

# Header
set title "Time evolution of the occupation of states x^2\n \
{Initalial state |2>, {/Symbol k} = 0.1E_0, \n \
Energies: E_0 = {/=12 @^{/Symbol=10 -}{/=14 h}}{/Symbol w}_0 = 1.0 meV \
and  {/=12 @^{/Symbol=10 -}{/=14 h}}{/Symbol W} = 0.4 meV"

# Key
set key right center 
set key font ",18"
set key spacing "1.8"

# Line colors
#      Yellow
set style line 2 lc rgb "#DB3502" lt 1 lw 1.5
set style line 3 lc rgb "#DB9A02" lt 1 lw 1.5
set style line 4 lc rgb "#EBE710" lt 1 lw 1.5
set style line 5 lc rgb "#F1F740" lt 1 lw 1.5
#      Green
set style line 6 lc rgb "#97E384" lt 1 lw 1.5
set style line 7 lc rgb "#9EE88B" lt 1 lw 1.5
set style line 8 lc rgb "#A5E895" lt 1 lw 1.5
set style line 9 lc rgb "#A5EB94" lt 1 lw 1.5

# Plot 
plot for [col=2:5] "Dataset1_kappa0.1/nit.dtx"  using 1:col with lines ls col title columnheader,\
#     for [col=2:5] "Dataset2_kappa0.05/nit.dtx" using 1:col with lines ls col   title columnheader,\
