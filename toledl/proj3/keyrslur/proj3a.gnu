reset

set terminal postscript portrait enhanced color solid defaultplex "Helvetica" 14
set output 'graphproj3a.eps'

#set xrange [0:256]
#set yrange [-35:0]
set xlabel "i"
set ylabel "log(abs(E_i)) both real and imag"
set title "Eigenspectrum"
set size 2,1
set key right center 
set key font ",18"
set key spacing "1.8"

plot "reigval.dtx" with lines lt rgb "#2C5A8F" title "Real part",\
     "ieigval.dtx"  lt rgb "#006400" title "Imaginary part",\
