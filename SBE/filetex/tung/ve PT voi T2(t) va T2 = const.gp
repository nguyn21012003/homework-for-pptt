# Đặt kích thước cửa sổ hiển thị và font chữ
set term wxt size 1920,1080 enhanced font 'Arial,12'

# Đặt nhãn cho các trục
set xlabel 'Time (fs)' font 'Arial,20'
set ylabel 'Total Polarization P(t)' font 'Arial,20'

# Chỉnh vị trí và font chữ cho chú thích
set key top right font "Arial,10"
set key font ",25"
set xtics font ",15"  
set ytics font ",15"
# Đặt phạm vi cho trục x và y
#set xrange [-100:200]   
#set yrange [0:1]     

# Hiển thị lưới
set grid
# Vẽ đồ thị 2D từ 5 file dữ liệu với các giá trị khác nhau của khi0
plot 'Nt&Pt_khi0_0.1_deltat_10_delta0_20_T2_50_time_1414.txt' using 1:3 with lines lw 2.5 lc rgb "red" title 'T_2 = 50', \
     'Nt&Pt_khi0_0.1_deltat_10_delta0_20_T2_100_time_1414.txt' using 1:3 with lines lw 2.5 lc  rgb "green" title 'T_2 = 100',\
     'Nt&Pt_khi0_0.1_deltat_10_delta0_20_T2_210_time_1415.txt' using 1:3 with lines  lw 2.5 lc rgb "orange" title 'T_2 = 210',\
     'Nt&Pt_khi0_0.1_deltat_10_delta0_20_T2(t)_time_1415.txt' using 1:3 with lines lw 2.5 lc rgb "blue" title 'T_2(t)'
     
