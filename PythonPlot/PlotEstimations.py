import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys #для приема имени файла для построения через аргумент команды

# data =pd.read_csv(sys.argv[1], delimiter=",");
data =pd.read_csv("data/Kalman_est.csv", delimiter=",");

time = np.linspace(0, (data.iloc[:, 0].size - 1)/100/60, data.iloc[:, 0].size) #в минутах /100/60

round = 100

'''Оценки Линейные скорости'''
fig1, (fig1_ax1, fig1_ax2) = plt.subplots(1, 2)
# Восточная составляющая
fig1.suptitle("Оценки скоростей")
fig1_ax1.set_title("Оценка восточной скорости");
fig1_ax1.set_xlabel("мин")
fig1_ax1.set_ylabel("м/с")
fig1_ax1.legend(loc="best")
fig1_ax1.plot(time, np.round(data["Ve"], round), label="Ve");
fig1_ax1.grid(True)
# Северная составляющая
fig1_ax2.set_title("Оценка северной скорости");
fig1_ax2.set_xlabel("мин")
fig1_ax2.set_ylabel("м/с")
fig1_ax2.legend(loc="best")
fig1_ax2.plot(time, np.round(data["Vn"],round), label="Vn");
fig1_ax2.grid(True)


'''Оценки углов ориентации'''
fig2, (fig2_ax1) = plt.subplots(1, 1)
fig2.suptitle("Оценки ориентации Фx, Фy")
# Угол Фx
fig2_ax1.set_xlabel("мин")
fig2_ax1.set_ylabel("угл. мин.")
fig2_ax1.plot(time, 60*np.round(np.rad2deg(data["Phi_e"]),round), label="Фx");
# Угол Фy
fig2_ax1.plot(time, 60*np.round(np.rad2deg(data["Phi_n"]),round), label="Фy");
fig2_ax1.legend(loc="best")
fig2_ax1.grid(True)

'''Осредняем скорости дрейфов гироскопов и сглаживаем ФНЧ с T=100'''

Tf = 1e7;
mean = [[0],[0]]
flf = [[0], [0]] # filter lower frequency
for i, dx in enumerate(data["d_omega_x"]):
    #MeanAb[i] = (Ldoub) MeanAb[i] + (Ab[i] - MeanAb[i]) / (iter + 1); //(iter * MeanAb[i] + Ab[i])/(iter + 1);
    mean[0].append(mean[0][i] + (dx - mean[0][i])/(i + 1));
    mean[1].append(mean[1][i] + (data["d_omega_y"][i+1] - mean[1][i])/(i + 1));
    flf[0].append((dx + Tf/100*flf[0][-1])/(1 + Tf/100))
    flf[1].append((data["d_omega_y"][i+1] + Tf/100*flf[1][-1])/(1 + Tf/100))


'''Оценки дрейфов гироскопов'''
fig3, (fig3_ax1) = plt.subplots(1, 1)
# Восточное направление
fig3.suptitle("Оценки скоростей дрейфов гироскопов")
fig3_ax1.set_xlabel("мин")
fig3_ax1.set_ylabel("град/час")
fig3_ax1.plot(time, np.round(np.rad2deg(data["d_omega_x"]), round)*3600, label="$\delta\omega_x$");
fig3_ax1.plot(time, np.round(np.rad2deg(data["d_omega_y"]), round)*3600, label="$\delta\omega_y$");
fig3_ax1.plot(time, np.round(np.rad2deg(mean[0][:-1]), round)*3600, label="$E(\omega_x)$");
fig3_ax1.plot(time, np.round(np.rad2deg(mean[1][:-1]), round)*3600, label="$E(\omega_y)$");
#ФНЧ
fig3_ax1.plot(time, np.round(np.rad2deg(flf[0][:-1]), round)*3600, label="$\omega_x^F$");
fig3_ax1.plot(time, np.round(np.rad2deg(flf[1][:-1]), round)*3600, label="$\omega_y^F$");
fig3_ax1.grid(True)
fig3_ax1.legend(loc="best")

plt.show()
