import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# data = pd.read_csv("/home/nikita_viksne/InertialNavigation/Data_files/TERM20250617091439765.log",
#                     delimiter=";", index_col = False, header = None);
# data.columns = ["Timestamp", "a_bx", "a_by", "a_bz",
#                 "Om_bx", "Om_by", "Om_bz"];

data = pd.read_csv("/home/nikita_viksne/InertialNavigation/Modelling_sensetive_elements/Data_files/data_acc_veloc_0_heading_50_freq_400_gps.csv",
                    delimiter=";", index_col = False, header = None);
data.columns = ["a_bx", "a_by", "a_bz",
                "Om_bx", "Om_by", "Om_bz",
                "BiasA_bx", "BiasA_by", "BiasA_bz", "BiasOm_bx", "BiasOm_by", "BiasOm_bz",
                "RandA_bx", "RandA_by", "RandA_bz", "RandOm_bx", "RandOm_by", "RandOm_bz",
                "Ve_gps", "Vn_gps", "RandVe_gps", "RandVn_gps", "bla-bla"];
data = data.drop("bla-bla", axis = 1)
time = np.linspace(0, (data.iloc[:, 0].size - 1)/400/60, data.iloc[:, 0].size) #в минутах /400/60 # (400 Hz)
meanOmb = [[0],[0],[0]]
meanAb = [[0],[0],[0]]

round = 100 #кооличество знаков после запятой

for i in range(3):# 3 оси
    for j in range(len(data["a_bx"])): # по времени
        # meanOmb[i].append(meanOmb[i][j] + (data[data.columns[3+i]][j] + data[data.columns[9+i]][j] + data[data.columns[15+i]][j] - meanOmb[i][j])/(j + 1));
        meanAb[i].append(meanAb[i][j] + (data[data.columns[i]][j] + data[data.columns[6+i]][j] + data[data.columns[12+i]][j] - meanAb[i][j])/(j + 1));
    del(meanOmb[i][0]) #удаляем первые начальные значения для сохранения длины массива
    del(meanAb[i][0]) #удаляем первые начальные значения для сохранения длины массива

# Акселерометры
fig5, (fig5_ax1, fig5_ax2, fig5_ax3) = plt.subplots(3, 1)
# Акселерометр X
fig5.suptitle("Гироскопы")
fig5_ax1.set_title("Показания акселерометра X");
fig5_ax1.set_xlabel("мин")
fig5_ax1.set_ylabel("м/(с^2)")
fig5_ax1.plot(time, np.round(data["a_bx"] + data["BiasA_bx"] + data["RandA_bx"], round), label="Показания");
fig5_ax1.plot(time, np.round(meanAb[0], round), linestyle='--', label="Среднее значение");
fig5_ax1.plot(time, np.round(data["a_bx"], round), label="Истинное занчение");
fig5_ax1.legend(loc="best")
fig5_ax1.grid(True)
# Акселерометр Y
fig5_ax2.set_title("Показания акселерометра Y");
fig5_ax2.set_xlabel("мин")
fig5_ax2.set_ylabel("м/(с^2)")
fig5_ax2.plot(time, np.round(data["a_by"] + data["BiasA_by"] + data["RandA_by"], round), label="Показания");
fig5_ax2.plot(time, np.round(meanAb[1], round), linestyle='--', label="Среднее значение");
fig5_ax2.plot(time, np.round(data["a_bx"], round), label="Истинное занчение");
fig5_ax2.legend(loc="best")
fig5_ax2.grid(True)
# Акселерометр Z
fig5_ax3.set_title("Показания акселерометра Z");
fig5_ax3.set_xlabel("мин")
fig5_ax3.set_ylabel("м/(с^2)")
fig5_ax3.plot(time, np.round(data["a_bz"] + data["BiasA_bz"] + data["RandA_bz"], round), label="Показания");
fig5_ax3.plot(time, np.round(meanAb[2], round), linestyle='--', label="Среднее значение");
fig5_ax3.plot(time, np.round(data["a_bz"], round), label="Истинное занчение");
fig5_ax3.legend(loc="best")
fig5_ax3.grid(True)

plt.show()


exit()

# Гироскопы
fig1, (fig1_ax1, fig1_ax2, fig1_ax3) = plt.subplots(3, 1)
# Гироскоп X
fig1.suptitle("Гироскопы")
fig1_ax1.set_title("Показания гироскопа X");
fig1_ax1.set_xlabel("мин")
fig1_ax1.set_ylabel("град/ч")
fig1_ax1.plot(time, np.round(np.rad2deg(data["Om_bx"] + data["BiasOm_bx"] + data["RandOm_bx"])*3600, round), label="Показания");
fig1_ax1.plot(time, np.round(np.rad2deg(data["Om_bx"])*3600, round), label="Истинное значение");
fig1_ax1.plot(time, np.round(np.rad2deg(meanOmb[0])*3600, round), linestyle='--', label="Среднее значение");

fig1_ax1.legend(loc="best")
fig1_ax1.grid(True)
# Гироскоп Y
fig1_ax2.set_title("Показания гироскопа Y");
fig1_ax2.set_xlabel("мин")
fig1_ax2.set_ylabel("град/ч")
fig1_ax2.plot(time, np.round(np.rad2deg(data["Om_by"] + data["BiasOm_by"] + data["RandOm_by"])*3600, round), label="Показания");
fig1_ax2.plot(time, np.round(np.rad2deg(data["Om_by"])*3600, round), label="Истинное значение");
fig1_ax2.plot(time, np.round(np.rad2deg(meanOmb[1])*3600, round), linestyle='--', label="Среднее значение");
fig1_ax2.legend(loc="best")
fig1_ax2.grid(True)
# Гироскоп Z
fig1_ax3.set_title("Показания гироскопа Z");
fig1_ax3.set_xlabel("мин")
fig1_ax3.set_ylabel("град/ч")
fig1_ax3.plot(time, np.round(np.rad2deg(data["Om_bz"] + data["BiasOm_by"] + data["RandOm_by"])*3600, round), label="Показания");
fig1_ax3.plot(time, np.round(np.rad2deg(data["Om_bz"])*3600, round), label="Истинное значение");
fig1_ax3.plot(time, np.round(np.rad2deg(meanOmb[2])*3600, round), linestyle='--', label="Среднее значение");
fig1_ax3.legend(loc="best")
fig1_ax3.grid(True)

plt.show()
