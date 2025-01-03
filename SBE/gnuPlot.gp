

    #set multiplot layout 1,2


    set ylabel "y"
    set xlabel "x"
    set zlabel "z"


    set grid
    set key horiz
    
    #splot "fe_ChiKhacNhau_deltaT=100_delta0=100at1010.data.sbe.txt" u 1:2:3 with lines tit "fe_t" 
    
    
    #set ylabel "y"
    #set xlabel "x"
    #set zlabel "z"
    #splot "N(t)_Chi=0.01_deltaT=100_delta0=100at1010.data.sbe.txt" u 1:2:4 with lines tit "p_t"


    set ylabel "E(V/m)" font "Arial, 15"
    set xlabel "Energy ħω (MeV)" font "Arial, 15"
    set key font ",25"
set xtics font ",15"
set ytics font ",15"
    set xrange [-50:50]
    #plot "N(t)_Chi=0.01_deltaT=100_delta0=100at1010.data.sbe.txt" u 2:3 with lines tit "Phân cực toàn phần"
    #plot "N(t)_Chi=0.01_deltaT=100_delta0=100at1010.data.sbe.txt" u 2:4 with lines tit "Mật độ toàn phần"
    plot "Alpha_ChiKhacNhau_deltaT=100_delta0=100at1010.data.sbe.txt" u 2:8 with lines lw 5 tit "E(ω)"
    #plot "N(t)_Chi=0.01_deltaT=100_delta0=100at1010.data.sbe.txt" u 1:3 with lines 



