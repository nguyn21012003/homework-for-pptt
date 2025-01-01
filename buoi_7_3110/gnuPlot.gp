
  set datafile separator "|"  
  splot "Laplace.txt" u 1:2:3 w l 
  unset multiplot
  pause -1

