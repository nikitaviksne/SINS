#!/bin/bash
# по новой версии программы, аргументы передаются в следующем порядке: курс в град, крен в град, тангаж в град, начальная скорость в м/с, файл откуда читать значения, файл, куда записывать значения навигационных параметров (опционально, если не переданно, будет выведено на печать). Курс крен и тангаж нужны для вычисления истинных значений навигационных параметров

build/debug/main 50 0 0 0 /home/nikita_viksne/Документы/Python/Modelling_sensetive_elements/Data_files/data_acc_veloc_0_heading_50_freq_400_gps.bin
#build/debug/main 50 0 0 30 /home/nikita_viksne/Документы/Python/Modelling_sensetive_elements/Data_files/data_acc_veloc_30_heading_50_freq_400_gps.bin
#build/debug/main 90 0 0 30 /home/nikita_viksne/Документы/Python/Modelling_sensetive_elements/Data_files/data_acc_veloc_30_heading_90_freq_400_gps.bin
