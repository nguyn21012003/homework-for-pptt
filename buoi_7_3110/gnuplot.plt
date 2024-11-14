set ylabel "y"
set xlabel "x"
set zlabel "V"

set grid
# set key horiz

splot "ElectricPotentials.txt" u 2:3:4 with lines tit "V"