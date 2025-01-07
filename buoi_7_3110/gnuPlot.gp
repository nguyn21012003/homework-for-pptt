
    set ylabel "y"
    set xlabel "x"
    set zlabel "V"


    set multiplot layout 1,2
    set grid
    # set key horiz

    splot "ElectricPotentials.txt" u 2:3:4 with lines tit "Jacobian"
    
    
    set ylabel "y"
    set xlabel "x"
    set zlabel "V"
    splot "ElectricPotentials.txt" u 2:3:5 with lines tit "Gaussian"


    unset multiplot

    pause -1
