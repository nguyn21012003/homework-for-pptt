
    set xlabel "x"
    set ylabel "y"
    set grid 
    set datafile separator ","

    plot "dataCircle.csv" u 3:4 with points ps 1 pt 7,    "dataCircle_out.csv" u 1:2 with points ps 1 pt 7
    pause -1

            