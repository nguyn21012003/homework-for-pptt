set multiplot layout 2, 1
set datafile separator ","

splot "HyperbolicEquationData.txt" using 1:2:3 with lines ,\
"HyperbolicEquationDataFric.txt" u 1:2:3 with lines

unset multiplot
pause -1