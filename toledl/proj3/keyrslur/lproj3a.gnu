reset

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'graphproj3alog.eps'



set xrange [0:256]
set yrange [-37:5]
set y2range [-37:0]
set ytics 5 nomirror tc lt 7  # lc rgb '#4CA305'
set y2tics 5 nomirror tc lt 6 # lc rgb '#A30D05'

set xzeroaxis

set xlabel "i"
set ylabel "log|Re(E_i)|"
set y2label "log|Im(E_i)|"
set title "Eigenvalues E_i in a system with damping strength {/Symbol k} = 0.2E_0 \
\n Initalial state |0>, Energies: E_0 = {/=12 @^{/Symbol=10 -}{/=14\
h}}{/Symbol w} = 1.0 meV and  {/=12 @^{/Symbol=10 -}{/=14 h}}{/Symbol W} = 0.4 meV"
set size 3,1
set key right center 
set key font ",18"
set key spacing "1.8"

set boxwidth 0.35
set style fill solid

plot "eigvallog.dtx" using 1:3 axes x1y2 with boxes lt 6 notitle,\
     "eigvallog.dtx" using 1:2 axes x1y1 lt 7 notitle,\
