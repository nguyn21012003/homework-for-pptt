

    #set multiplot layout 1,2


    set ylabel "y"
    set xlabel "x"
    set zlabel "z"
    set xrange [-100:100]
    

    set grid
    set key horiz

    #splot "25&10fs with N=100&0.5.data.sbe.txt" u 1:2:3 with lines tit "fe_t" 
    
    
    #set ylabel "y"
    #set xlabel "x"
    #set zlabel "z"
    #splot "25&10fs with N=100&0.5.data.sbe.txt" u 1:2:10 with lines tit "p_t"


    set xlabel "Omega"
    set ylabel "alpha(omega)"
    #plot "25&10fs with N=100&0.5.data.absortion.txt" u 2:3 with lines tit "Phân cực toàn phần"
    #plot "25&10fs with N=100&0.5.data.absortion.txt" u 2:10 with lines tit "Mật độ toàn phần"
    plot "Alpha_ChiKhacNhau_deltaT=20_delta0=10at1632.data.sbe.txt" u 2:5 with lines tit "chi = 0.001" ,\
         "Alpha_ChiKhacNhau_deltaT=20_delta0=10at1633.data.sbe.txt" u 2:5 with lines tit "chi = 0.01" ,\
         "Alpha_ChiKhacNhau_deltaT=10_delta0=20at1627.data.sbe.txt" u 2:5 with lines tit "chi = 0.1" ,\
         "Alpha_ChiKhacNhau_deltaT=10_delta0=20at1628.data.sbe.txt" u 2:5 with lines tit "chi = 0.2" ,\
         "Alpha_ChiKhacNhau_deltaT=10_delta0=20at1629.data.sbe.txt" u 2:5 with lines tit "chi = 0.5" ,\
         "Alpha_ChiKhacNhau_deltaT=20_delta0=10at1634.data.sbe.txt" u 2:5 with lines tit "chi = 2" ,\



    #unset multiplot

    pause -1

