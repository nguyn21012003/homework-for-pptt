

    #set multiplot layout 1,3


    set ylabel "y"
    set xlabel "x"
    set zlabel "z"


    set grid
    #set key horiz

    splot "heatEquationData.txt" u 1:2:3 with lines 

    set datafile separator '	'

    #set ylabel "y"
    #set xlabel "x"
    #set zlabel "z"
    #splot "heatEquationData.txt" u 1:2:4 with lines tit "p_t"


    #set xlabel "Omega"
    #set ylabel "alpha(omega)"
    #plot "heatEquationData.txt" u 2:3 with lines tit "Phân cực toàn phần"
    ##plot "heatEquationData.txt" u 2:4 with lines tit "Mật độ toàn phần"
    #plot "heatEquationData.txt" u 2:5 with lines tit "Phổ hấp thụ"


    #unset multiplot

    pause -1

