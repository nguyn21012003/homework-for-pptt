# Đặt kích thước cửa sổ hiển thị và font chữ
set term wxt size 1920,1080 enhanced font 'Arial,12'

# Đặt nhãn cho các trục
set xlabel 'Năng lượng (meV)' font 'Arial,14'
set ylabel 'Phổ hấp thụ' font 'Arial,20'

# Chỉnh vị trí và font chữ cho chú thích
set key top right font "Arial,12"
set key font ",25"
set xtics font ",15"  
set ytics font ",15"
# Đặt phạm vi cho trục x và y
set xrange [-100:50]   
#set yrange [-500:500]     

# Hiển thị lưới
set grid

# Vẽ đồ thị 2D từ 5 file dữ liệu với các giá trị khác nhau của khi0
plot 'alpha(omega)_khi0_0.001_deltat_10_delta0_20_T2_210_time_1425.txt' using 1:4 with lines lw 3 lc rgb "red" title 'χ_0 = 0.001', \
     'alpha(omega)_khi0_0.1_deltat_10_delta0_20_T2_210_time_1428.txt' using 1:4 with lines lw 3 lc rgb "green" title 'χ_0 = 0.05', \
     'alpha(omega)_khi0_0.2_deltat_10_delta0_20_T2_210_time_1431.txt' using 1:4 with lines lw 3 lc rgb "orange" title 'χ_0 = 0.2',\
     'alpha(omega)_khi0_0.5_deltat_10_delta0_20_T2_210_time_1435.txt' using 1:4 with lines lw 3 lc rgb "blue" title 'χ_0 = 0.5'
