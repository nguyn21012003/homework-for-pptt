
    set xlabel "x"
    set ylabel "y"
    set zlabel "z"


    set grid
    set key horiz

    plot "Jacobi&GaussSeidel.txt" u 1:2 with lines 



    unset multiplot

    pause -1

