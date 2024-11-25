

    #set multiplot layout 1,3


    set ylabel "y"
    set xlabel "x"
    set zlabel "z"


    set grid
    #set key horiz

    splot "heatEquationData.txt" u 1:2:3 with lines 

    #set datafile separator '	'
    

    #unset multiplot

    pause -1

