import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys #для приема имени файла для построения через аргумент команды

path = sys.argv[1]
data =pd.read_csv("~/InertialNavigation/"+path, delimiter=",");

time = np.linspace(0, (data.iloc[:, 0].size - 1)/100/60, data.iloc[:, 0].size) #в минутах /100/60, где 100 это частота выдачи оценок

round = 100

'''Оценки линейных скоростей'''
fig1, (fig1_ax1, fig1_ax2) = plt.subplots(1, 2)
# Восточная составляющая
fig1.suptitle("Скорости")
fig1_ax1.set_title("Оценка восточной скорости");
fig1_ax1.set_xlabel("мин")
fig1_ax1.set_ylabel("м/с")
fig1_ax1.plot(time, np.round(data["Ve"], round), label="$\hat{V}e$");
fig1_ax1.legend(loc="best")
fig1_ax1.grid(True)
# Северная составляющая
fig1_ax2.set_title("Оценка северной скорости");
fig1_ax2.set_xlabel("мин")
fig1_ax2.set_ylabel("м/с")
fig1_ax2.plot(time, np.round(data["Vn"],round), label="$\hat{V}n$");
fig1_ax2.legend(loc="best")
fig1_ax2.grid(True)

'''Углы платформы (ориентация)'''
fig2, (fig2_ax1) = plt.subplots(1, 1)
fig2.suptitle("Оценки углов платформы")
# Уголы
fig2_ax1.set_xlabel("сек")
fig2_ax1.set_ylabel("град")
fig2_ax1.plot(time, np.round(3600*np.rad2deg(data["Phi_e"]),round), label="$\hat{\Phi}_x$");
fig2_ax1.plot(time, np.round(3600*np.rad2deg(data["Phi_n"]),round), label="$\hat{\Phi}_n$");
fig2_ax1.legend(loc="best")
fig2_ax1.grid(True)


'''Оценки скоростей дрейфов'''
fig3, (fig3_ax1) = plt.subplots(1, 1)
# Восточное направление
fig3.suptitle("Оценка скоростей дрейфов гироскопов")
fig3_ax1.set_xlabel("мин")
fig3_ax1.set_ylabel("град/ч")
fig3_ax1.plot(time, np.round(data["d_omega_x"], round), label="$\hat{\delta\omega}_x$");
fig3_ax1.plot(time, np.round(data["d_omega_y"], round), label="$\hat{\delta\omega}_y$");
fig3_ax1.legend(loc="best")
fig3_ax1.grid(True)

plt.show()
