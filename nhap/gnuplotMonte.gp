set datafile separator ","

plot "circleIN.csv" u 1:2 with points ps 1 pt 7,\
    "circleOUT.csv" u 1:2 with points ps 1 pt 7