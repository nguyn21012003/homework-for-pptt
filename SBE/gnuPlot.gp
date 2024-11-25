

    set multiplot layout 1,2


    set ylabel "y"
    set xlabel "x"
    set zlabel "z"


    set grid
    set key horiz

    #splot "50&50fs with N=100&0.5.data.sbe.txt" u 1:2:3 with lines tit "fe_t"
    
    
    set ylabel "y"
    set xlabel "x"
    set zlabel "z"
    #splot "50&50fs with N=100&0.5.data.sbe.txt" u 1:2:4 with lines tit "p_t"


    set xlabel "Omega"
    set ylabel "alpha(omega)"
    plot "50&50fs with N=100&0.5.data.absortion.txt" u 2:3 with lines tit "Phân cực toàn phần"
    plot "50&50fs with N=100&0.5.data.absortion.txt" u 2:4 with lines tit "Mật độ toàn phần"
    #plot "50&50fs with N=100&0.5.data.absortion.txt" u 2:5 with lines tit "Phổ hấp thụ"


    unset multiplot

    pause -1

