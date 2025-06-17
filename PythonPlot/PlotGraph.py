import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys #для приема имени файла для построения через аргумент команды

relative_path = sys.argv[1]
# data =pd.read_csv("~/InertialNavigation/data/"+relative_path, delimiter=";");
data =pd.read_csv(relative_path, delimiter=";");


time = np.linspace(0, (data.iloc[:, 0].size - 1)/25/60, data.iloc[:, 0].size) #в минутах /100/60

round = 100

'''Линейные скорости'''
fig1, (fig1_ax1, fig1_ax2) = plt.subplots(1, 2)
# Восточная составляющая
fig1.suptitle("Скорости")
fig1_ax1.set_title("Восточная скорость");
fig1_ax1.set_xlabel("мин")
fig1_ax1.set_ylabel("м/с")
fig1_ax1.legend(loc="best")
fig1_ax1.plot(time, np.round(data["Ve"], round), label="Ve");
fig1_ax1.grid(True)
# Северная составляющая
fig1_ax2.set_title("Северная скорость");
fig1_ax2.set_xlabel("мин")
fig1_ax2.set_ylabel("м/с")
fig1_ax2.legend(loc="best")
fig1_ax2.plot(time, np.round(data["Vn"],round), label="Vn");
fig1_ax2.grid(True)

'''Углы ориентации'''
fig2, (fig2_ax1, fig2_ax2, fig2_ax3) = plt.subplots(3, 1)
fig2.suptitle("Углы")
# Угол курса
fig2_ax1.set_title("Угол курса");
fig2_ax1.set_xlabel("мин")
fig2_ax1.set_ylabel("угл. мин")
fig2_ax1.legend(loc="best")
fig2_ax1.plot(time, np.round(np.rad2deg(data["Heading"])*60,round), label="Курс");
fig2_ax1.grid(True)
# Угол крена
fig2_ax2.set_title("Угол крена");
fig2_ax2.set_xlabel("мин")
fig2_ax2.set_ylabel("угл. мин")
fig2_ax2.legend(loc="best")
fig2_ax2.plot(time, np.round(np.rad2deg(data["Roll"])*60,round), label="Крен");
fig2_ax2.grid(True)
# Угол тангажа
fig2_ax3.set_title("Угол тангажа");
fig2_ax3.set_xlabel("мин")
fig2_ax3.set_ylabel("угл. мин")
fig2_ax3.legend(loc="best")
fig2_ax3.plot(time, np.round(np.rad2deg(data["Pitch"])*60,round), label="Тангаж");
fig2_ax3.grid(True)

'''Ошибки по координатам'''
fig3, (fig3_ax1, fig3_ax2) = plt.subplots(1, 2)
# Восточное направление
fig3.suptitle("Ошибка по координатам")
fig3_ax1.set_title("Восточное направление");
fig3_ax1.set_xlabel("мин")
fig3_ax1.set_ylabel("м")
fig3_ax1.legend(loc="best")
fig3_ax1.plot(time, np.round(data["d_E"], round), label="dE");
fig3_ax1.grid(True)
# Северное направление
fig3_ax2.set_title("Северное направление");
fig3_ax2.set_xlabel("мин")
fig3_ax2.set_ylabel("м")
fig3_ax2.legend(loc="best")
fig3_ax2.plot(time, np.round(data["d_N"],round), label="dN");
fig3_ax2.grid(True)

'''Траектория'''
fig3, (fig3_ax1) = plt.subplots(1, 1)
# Восточное направление
fig3.suptitle("Траектория")
fig3_ax1.set_xlabel("град")
fig3_ax1.set_ylabel("град")
fig3_ax1.plot(np.round(np.rad2deg(data["Lambda"]), round), np.round(np.rad2deg(data["Phi"]), round) );
fig3_ax1.grid(True)

plt.show()
