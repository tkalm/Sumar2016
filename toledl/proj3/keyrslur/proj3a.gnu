reset

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'graphproj3a.eps'

set xrange [0:256]
set yrange [-20:20]
set y2range [-1.4:0]
set ytics 5 nomirror tc lt 1
set y2tics 0.2 nomirror tc lt 2

#set logscale y
set xlabel "i"
set ylabel "Re(E_i)"
set y2label "Im(E_i)"
set title "Eigenspectrum"
set size 2,1
set key right center 
set key font ",18"
set key spacing "1.8"

plot "eigval.dtx" using 1:3 axes x1y2 with boxes lt 2 notitle,\
     "eigval.dtx" using 1:2 axes x1y1 lt 1 notitle,\
