

    #set multiplot layout 1,3

    set datafile separator '	'

    set title "Heat Equation Solution"
    set xlabel "xnguyendz"
    set ylabel "ynguyendz"
    set zlabel "T(x,y)"

    set grid
    #set key horiz

    splot "heatEquationData.txt" u 1:2:3 with lines 

    #set datafile separator '	'
    

    #unset multiplot

    pause -1

