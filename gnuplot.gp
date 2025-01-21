

    #set multiplot layout 1,2


    set ylabel "y"
    set xlabel "x"
    set zlabel "z"


    set grid
    set key horiz

    #splot "fe_ChiKhacNhau_deltaT=100_delta0=100at1552.data.sbe.txt" u 1:2:3 with lines tit "fe_t" 
    
    
    #set ylabel "y"
    #set xlabel "x"
    #set zlabel "z"
    #splot "N(t)_Chi=2_deltaT=100_delta0=100at1552.data.sbe.txt" u 1:2:4 with lines tit "p_t"


    set xlabel "Omega"
    set ylabel "alpha(omega)"
    #plot "N(t)_Chi=2_deltaT=100_delta0=100at1552.data.sbe.txt" u 2:3 with lines tit "Phân cực toàn phần"
    #plot "N(t)_Chi=2_deltaT=100_delta0=100at1552.data.sbe.txt" u 2:4 with lines tit "Mật độ toàn phần"
    #plot "Alpha_ChiKhacNhau_deltaT=100_delta0=100at1552.data.sbe.txt" u 2:5 with lines 
    plot "N(t)_Chi=2_deltaT=100_delta0=100at1552.data.sbe.txt" u 1:3 with lines 


    #unset multiplot

    pause -1

