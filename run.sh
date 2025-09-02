#!/bin/bash
# по новой версии программы, аргументы передаются в следующем порядке: курс в град, крен в град, тангаж в град, начальная скорость в м/с, файл откуда читать значения, постоянные погрешности акселерометров (1 включены, 0 выключены), постоянные погрешности гироскопов, случайные погрешности акселерометров, случайные погрешности гироскопов, случайные погрешности корректора (по скорости). Курс крен и тангаж нужны для вычисления истинных значений навигационных параметров

#build/debug/main 50 0 0 0 /home/nikita/Документы/Python/Modelling_sensetive_elements/Data_files/data_acc_veloc_0_heading_50_freq_400_gps.csv
#build/debug/main 50 0 0 0 /home/nikita/Документы/Python/Modelling_sensetive_elements/Data_files/data_acc_veloc_0_heading_50_freq_400_gps.bin 1 1 1 1 1 
#build/debug/main 50 0 0 30 /home/nikita/Документы/Python/Modelling_sensetive_elements/Data_files/data_acc_veloc_30_heading_50_freq_400_gps.bin
build/debug/main 90 0 0 30 /home/nikita/Документы/C_Cpp_progs/InertialNavigation/Data_files/data_acc_veloc_30_heading_90_freq_100_V_coo_gps.csv 1 1 1 1 1 1

# Поворот
#build/main 90 0 0 30 /home/nikita/Документы/C_Cpp_progs/InertialNavigation/Data_files/data_acc_veloc_30_heading_90_freq_400_turn_V_coo_gps.csv 1 1 1 1 1 1

#Реальные сырые данные
# build/debug/main 50 0 0 0 /home/nikita/Документы/C_Cpp_progs/InertialNavigation/Data_files/TERM20250616144529640.log 1 1 1 1 1 
