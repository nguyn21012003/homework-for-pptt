# Gnuplot script to plot Gaussian pulse
set xrange[-250 to 250]
set xlabel "Time (fs)" font "Arial, 15" offset 0,-1
set ylabel "E(V/m)" font "Arial, 15" 
set key font ",25"
set xtics font ",15"
set ytics font ",15"
set grid
plot "gaussian_pulse_data.txt" u 1:2 w l lw 5 tit "E(V/m)"

