set datafile separator ","

set terminal pngcairo size 800,600
set terminal pdf
set output 'numericalmethod.gnuplot.pdf'

set title r"Nghiệm của phương trình vi phân y'=y-t^2+1"
set xlabel "t"
set ylabel "y"
set key left top

# Plot the different methods
plot "data.txt" using 2:3 with lines title "Exact", \
     "data.txt" using 2:4 with linespoints title "Euler", \
     "data.txt" using 2:5 with linespoints title "Euler Modified", \
     "data.txt" using 2:6 with linespoints title "RK2", \
     "data.txt" using 2:7 with linespoints title "RK3", \
     "data.txt" using 2:8 with linespoints title "RK4"
