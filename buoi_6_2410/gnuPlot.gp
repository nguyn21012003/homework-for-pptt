
    set ylabel "y"
    set xlabel "x"
    set zlabel "z"


    set grid
    set key horiz

    set datafile separator ","
    

    plot "nonlinearOSC.txt" u 1:2 with lines tit "không có trường ngoài và ma sát",        "nonlinearOSC.txt" u 1:3 with lines tit "có trường ngoài và không có ma sát",        "nonlinearOSC.txt" u 1:4 with lines tit "không có trường ngoài và có ma sát",        "nonlinearOSC.txt" u 1:5 with lines tit "có trường ngoài và ma sát"


    unset multiplot

    pause -1

