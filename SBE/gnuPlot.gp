

    #set multiplot layout 1,2


    set ylabel "y"
    set xlabel "x"
    set zlabel "z"


    set grid
    set key horiz

    #splot "20&10fs with N=100&2.data.sbe.txt" u 1:2:3 with lines tit "fe_t" 
    
    
    #set ylabel "y"
    #set xlabel "x"
    #set zlabel "z"
    #splot "20&10fs with N=100&2.data.sbe.txt" u 1:2:4 with lines tit "p_t"


    set xlabel "Omega"
    set ylabel "alpha(omega)"
    #plot "20&10fs with N=100&2.data.absortion.txt" u 2:3 with lines tit "Phân cực toàn phần"
    #plot "20&10fs with N=100&2.data.absortion.txt" u 2:4 with lines tit "Mật độ toàn phần"
    plot "20&10fs with N=100&2.data.absortionT2change.txt" u 1:3 with lines  ,\
         "20&10fs with N=100&2.data.absortion.txt" u 1:3 with lines 

    


    #unset multiplot

    pause -1

