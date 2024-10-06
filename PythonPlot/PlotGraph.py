import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data =pd.read_csv("~/InertialNavigation/data/Nav_res.csv", delimiter=";");

time = np.linspace(0, data.iloc[:, 0].size - 1, data.iloc[:, 0].size) #в минутах /100/60

'''Линейные скорости'''
fig1, (fig1_ax1, fig1_ax2) = plt.subplots(1, 2)
# Восточная составляющая
fig1_ax1.set_title("Восточная скорость");
fig1_ax1.set_xlabel("мин")
fig1_ax1.set_ylabel("м/с")
fig1_ax1.legend(loc="best")
fig1_ax1.plot(time, np.round(data["Ve"], 1), label="Ve");
fig1_ax1.grid(True)
# Северная составляющая
fig1_ax2.set_title("Северная скорость");
fig1_ax2.set_xlabel("мин")
fig1_ax2.set_ylabel("м/с")
fig1_ax2.legend(loc="best")
fig1_ax2.plot(time, np.round(data["Vn"],1), label="Vn");
fig1_ax2.grid(True)

'''Углы ориентации'''
fig2, (fig2_ax1, fig2_ax2, fig2_ax3) = plt.subplots(3, 1)
# Угол курса
fig2_ax1.set_title("Угол курса");
fig2_ax1.set_xlabel("мин")
fig2_ax1.set_ylabel("град")
fig2_ax1.legend(loc="best")
fig2_ax1.plot(time, np.round(np.rad2deg(data["Heading"]),1), label="Курс");
fig2_ax1.grid(True)
# Угол крена
fig2_ax2.set_title("Угол крена");
fig2_ax2.set_xlabel("мин")
fig2_ax2.set_ylabel("град")
fig2_ax2.legend(loc="best")
fig2_ax2.plot(time, np.round(np.rad2deg(data["Roll"]),1), label="Крен");
fig2_ax2.grid(True)
# Угол тангажа
fig2_ax3.set_title("Угол тангажа");
fig2_ax3.set_xlabel("мин")
fig2_ax3.set_ylabel("град")
fig2_ax3.legend(loc="best")
fig2_ax3.plot(time, np.round(np.rad2deg(data["Pitch"]),1), label="Тангаж");
fig2_ax3.grid(True)

plt.show()