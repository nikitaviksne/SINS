import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys #для приема имени файла для построения через аргумент команды

path = sys.argv[1]
data = pd.read_csv(path, delimiter=";");

data.columns = ["Abx", "Aby", "Abz", "Ombx", "Omby", "Ombz", "biasAbx", "biasAby", "biasAbz", "biasOmbx", "biasOmby", "biasOmbz", "randAbx", "randAby", "randAbz", "randOmbx", "randOmby", "randOmbz", "Vgps_x", "Vgps_y", "randVgps_x", "randVgps_y", "Lat_gps", "Lon_gps", "randLat_gps", "randLon_gps", "bla-bla"]

del(data["bla-bla"])

time = np.linspace(0, (data.iloc[:, 0].size - 1)/400/60, data.iloc[:, 0].size) #в минутах /400/60
round = 100

'''Линейные скорости'''
fig1, (fig1_ax1, fig1_ax2) = plt.subplots(1, 2)
# Восточная составляющая
fig1.suptitle("Скорости")
fig1_ax1.set_title("Восточная составляющая");
fig1_ax1.set_xlabel("мин")
fig1_ax1.set_ylabel("м/с")
fig1_ax1.plot(time, np.round(data["Vgps_x"], round), label="Ve");
# fig1_ax1.plot(time, np.round(estimations["Ve"], round), linestyle='--', label="$\hat{Ve}$");
fig1_ax1.legend(loc="best")
fig1_ax1.grid(True)
# Северная составляющая
fig1_ax2.set_title("Северная составляющая");
fig1_ax2.set_xlabel("мин")
fig1_ax2.set_ylabel("м/с")
fig1_ax2.plot(time, np.round(data["Vgps_y"],round), label="Vn");
# fig1_ax2.plot(time, np.round(estimations["Vn"],round), linestyle='--', label="$\hat{Vn}$");
fig1_ax2.legend(loc="best")
fig1_ax2.grid(True)

'''Траектория'''
fig2, (fig2_ax1) = plt.subplots(1, 1)
fig2.suptitle("Траектория")
fig2_ax1.set_xlabel("град")
fig2_ax1.set_ylabel("град")
fig2_ax1.plot(np.round(np.rad2deg(data["Lon_gps"]), round),np.round(np.rad2deg(data["Lat_gps"]),round));
# fig1_ax1.plot(time, np.round(estimations["Ve"], round), linestyle='--', label="$\hat{Ve}$");
fig2_ax1.legend(loc="best")
fig2_ax1.grid(True)


'''Угловые скорости (показания с гироскопов)'''
fig3, (fig3_ax1, fig3_ax2) = plt.subplots(2,1)
fig3.suptitle("Угловые скорости в связанной С.К.")
fig3_ax1.set_xlabel("мин")
fig3_ax1.set_ylabel("град/ч");
fig3_ax1.plot(time, 3600*np.round(np.rad2deg(data["Ombx"]), round),label="Ombx")
fig3_ax1.plot(time, 3600*np.round(np.rad2deg(data["Omby"]), round),label="Omby")
fig3_ax1.legend(loc="best")
fig3_ax1.grid(True)
#
fig3_ax2.set_xlabel("мин")
fig3_ax2.set_ylabel("град/ч");
fig3_ax2.plot(time, 3600*np.round(np.rad2deg(data["Ombz"]), round),label="Ombz")
fig3_ax2.legend(loc="best")
fig3_ax2.grid(True)

'''Ускорения в связанной системе координат (показания акселерометров)'''
fig4, (fig4_ax1, fig4_ax2) = plt.subplots(2,1)
fig4.suptitle("Линейные ускорения в связанной С.К.")
fig4_ax1.set_xlabel("мин")
fig4_ax1.set_ylabel("м/с^2");
fig4_ax1.plot(time, np.round(data["Abx"], round),label="Abx")
fig4_ax1.plot(time, np.round(data["Aby"], round),label="Aby")
fig4_ax1.legend(loc="best")
fig4_ax1.grid(True)
#
fig4_ax2.set_xlabel("мин")
fig4_ax2.set_ylabel("м/с^2");
fig4_ax2.plot(time, np.round(data["Abz"], round),label="Abz")
fig4_ax2.grid(True)

plt.show()
