# Đặt kích thước cửa sổ hiển thị và font chữ
set term wxt size 1920,1080 enhanced font 'Arial,14'

# Đặt nhãn cho các trục
set xlabel 'Time (fs)' font 'Arial,20'
set ylabel 'Total Density N(t)' font 'Arial,20'

# Chỉnh vị trí và font chữ cho chú thích
set key top right font "Arial,12"
set key font ",25"
# Đặt phạm vi cho trục x và y
set xrange [-50:50]   
#set yrange [0:1]     

set xtics font ",15"  
set ytics font ",15"
set key at -30,5.9e18
# Hiển thị lưới
set grid

# Vẽ đồ thị 2D từ 5 file dữ liệu với các giá trị khác nhau của khi0
plot 'Nt_khi0_0.1_deltat_10_delta0_20_T2(t)_time_1353.txt' using 1:2 with lines lw 5 lc rgb "red" title 'χ_0 = 0.1', \
     'Nt_khi0_0.5_deltat_10_delta0_20_T2(t)_time_1359.txt' using 1:2 with lines lw 5 lc rgb "blue" title 'χ_0 = 0.5', \
     'Nt_khi0_2_deltat_10_delta0_20_T2(t)_time_1405.txt' using 1:2 with lines lw 5 lc rgb "green" title 'χ_0 = 2', \
     'Nt_khi0_3_deltat_10_delta0_20_T2(t)_time_1411.txt' using 1:2 with lines lw 5 lc rgb "orange" title 'χ_0 = 3'

