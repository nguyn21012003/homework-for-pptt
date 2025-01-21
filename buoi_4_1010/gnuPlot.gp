
    set ylabel "y"
    set xlabel "x"
    set zlabel "z"

    set datafile separator "|"

    set grid
    set key horiz

    plot "data.txt" u 1:2 with lines ,    "data.txt" u 1:3 with lines title "σ"

    
    pause -1

