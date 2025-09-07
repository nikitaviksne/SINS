import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

path = sys.argv[1] #Data_files/data_acc_veloc_30_heading_50_R_0_P_0_freq_400_turn_V_coo_gps.csv

data = pd.read_csv(path,  delimiter=";", index_col = False, header = None)

data.columns = ["Abx", "Aby", "Abz", "Ombx", "Omby", "Ombz", "biasAbx", "biasAby", "biasAbz", "biasOmbx", "biasOmby", "biasOmbz", "randAbx", "randAby", "randAbz", "randOmbx", "randOmby", "randOmbz", "Vgps_x", "Vgps_y", "randVgps_x", "randVgps_y", "phi_gps", "lambda_gps", "randPhi_gps", "randLambda_gps", "bla-bla"]

data.drop("bla-bla", axis = 1)

# нахождение частоты данных из имени файла
splitted = path.split("_") # разюиваем строку по разделителям
pos_freq = splitted.index("freq")# ищем  номер позиции с ключевым словом "freq"
freq = int(splitted[pos_freq + 1]) # следующий за pos_freq это и есть частота
time = np.linspace(0, data.shape[0] / freq / 60, data.shape[0] )

fig1, (fig1_ax1) = plt.subplots(ncols=1, nrows=1)
fig1.suptitle(f"Истинные скорости (GPS без погрешностей), частота данных {freq} Гц")
# fig1_ax1.set_title("");
fig1_ax1.plot(time, data["Vgps_x"])
fig1_ax1.plot(time, data["Vgps_y"])
fig1_ax1.grid(True)

# fig3, (fig3_ax1) = plt.subplots(ncols = 1, nrows = 1)
# fig3.suptitle("Углы ориентации (Эйлера)")
# fig3_ax1.plot(np.rad2deg(data["heading"]))

fig5, (fig5_ax1, fig5_ax2, fig5_ax3) = plt.subplots(ncols=1, nrows=3)
fig5.suptitle(f"Показания ДУС, частота данных {freq} Гц")
fig5_ax1.set_title("$\Omega_{bx}$")
fig5_ax1.plot(time, np.rad2deg(data["Ombx"])*3600)
fig5_ax1.grid(True)
# 
fig5_ax2.set_title("$\Omega_{by}$")
fig5_ax2.plot(time, np.rad2deg(data["Omby"])*3600)
fig5_ax2.grid(True)
# 
fig5_ax3.set_title("$\Omega_{bz}$")
fig5_ax3.plot(time, np.rad2deg(data["Ombz"])*3600)
fig5_ax3.grid(True)

'''Ускорения (проверка кориолиса в момент разворота)'''

fig10, (fig10_ax1, fig10_ax2, fig10_ax3) = plt.subplots(ncols=1, nrows=3)
fig10.suptitle(f"Показания акселерометров, частота данных {freq} Гц")
fig10_ax1.set_title("$A_{bx}$")
fig10_ax1.plot(time, data["Abx"])
fig10_ax1.grid(True)
# 
fig10_ax2.set_title("$A_{by}$")
fig10_ax2.plot(time, data["Aby"])
fig10_ax2.grid(True)
# 
fig10_ax3.set_title("$A_{bz}$")
fig10_ax3.plot(time, data["Abz"])
fig10_ax3.grid(True)

plt.show()