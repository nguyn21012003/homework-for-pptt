

    set multiplot layout 1,3


    set ylabel "y"
    set xlabel "x"
    set zlabel "z"


    set grid
    set key horiz

    splot "50&50fs with N=100&0.0001.data.sbe.txt" u 1:2:3 with lines tit "fe_t"
    
    
    set ylabel "y"
    set xlabel "x"
    set zlabel "z"
    splot "50&50fs with N=100&0.0001.data.sbe.txt" u 1:2:4 with lines tit "p_t"


    set xlabel "Omega"
    set ylabel "alpha(omega)"
    plot "50&50fs with N=100&0.0001.data.absortion.txt" u 2:3 with lines tit "Phổ hấp thụ"


    unset multiplot

    pause -1

