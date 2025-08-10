import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys #для приема имени файла для построения через аргумент команды

# data =pd.read_csv(sys.argv[1], delimiter=",");
data = pd.read_csv("data/Kalman_est.csv", delimiter=",");

time = np.linspace(0, (data.iloc[:, 0].size - 1)/100/60, data.iloc[:, 0].size) #в минутах /100/60

round = 100

'''Оценки ошибок линейный скоростей'''
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


'''Оценки ориентации'''
fig2, (fig2_ax1, fig2_ax2) = plt.subplots(1, 2)
fig2.suptitle("Оценки углов");
fig2_ax1.set_title("Оценки ориентации Фx, Фy")
# Угол Фx
fig2_ax1.set_xlabel("мин")
fig2_ax1.set_ylabel("угл. мин.")
fig2_ax1.plot(time, 60*np.round(np.rad2deg(data["Phi_e"]),round), label="Фx");
# Угол Фy
fig2_ax1.plot(time, 60*np.round(np.rad2deg(data["Phi_n"]),round), label="Фy");
fig2_ax1.legend(loc="best")
fig2_ax1.grid(True)
#
fig2_ax2.set_title("Оценки углов ориентации крен, тангаж")
# Угол gamma
fig2_ax2.set_xlabel("мин")
fig2_ax2.set_ylabel("угл. мин.")
fig2_ax2.plot(time, 60*np.round(np.rad2deg(data["d_Roll"]),round), label="$\Delta\gamma$");
# Угол theta
fig2_ax2.plot(time, 60*np.round(np.rad2deg(data["d_Pitch"]),round), label="$\Delta\\theta$");
fig2_ax2.legend(loc="best")
fig2_ax2.grid(True)

'''Осредняем скорости дрейфов гироскопов и сглаживаем ФНЧ с T=100'''

Tf = 1e6;
alpha_f = 1e-4;
mean = [[0],[0]] # среднее значение ФНЧ
flf = [[0], [0]] # filter lower frequency
flf_alpha = [[0], [0]] # filter lower frequency с альфой
for i, dx in enumerate(data["d_omega_x"]):
    #MeanAb[i] = (Ldoub) MeanAb[i] + (Ab[i] - MeanAb[i]) / (iter + 1); //(iter * MeanAb[i] + Ab[i])/(iter + 1);
    flf[0].append((dx + Tf/100*flf[0][-1])/(1 + Tf/100))
    flf[1].append((data["d_omega_y"][i+1] + Tf/100*flf[1][-1])/(1 + Tf/100))
    # Фильтр с альфой
    flf_alpha[0].append(flf_alpha[0][-1]*(1-alpha_f) + dx*alpha_f);
    flf_alpha[1].append(flf_alpha[1][-1]*(1-alpha_f) + data["d_omega_y"][i+1]*alpha_f);
    mean[0].append(mean[0][i] + (flf[0][i+1] - mean[0][i])/(i + 1));
    mean[1].append(mean[1][i] + (flf[1][i+1] - mean[1][i])/(i + 1));
del(mean[0][0])
del(mean[1][0])
del(flf[0][0])
del(flf[1][0])
del(flf_alpha[0][0])
del(flf_alpha[1][0])

'''Оценки дрейфов гироскопов'''
fig4, (fig4_ax1, fig4_ax2) = plt.subplots(2, 1)
# Восточное направление
fig4.suptitle("Оценка скоростей дрейфов гироскопов")
fig4_ax1.set_title("Оценка скоростей дрейфа гироскопа X")
fig4_ax1.set_xlabel("мин")
fig4_ax1.set_ylabel("град/час")
fig4_ax1.plot(time, np.round(np.rad2deg(data["d_omega_x"]), round)*3600, label="$\delta\omega_x$");
# fig4_ax1.plot(time, np.round(np.rad2deg(mean[0]), round)*3600, label="$E(\omega_x^F)$");
# fig4_ax1.plot(time, np.round(np.rad2deg(flf[0]), round)*3600, label="$\omega_x^F$", alpha = 0.5, linestyle = "--", marker = "*");
fig4_ax1.plot(time, np.round(np.rad2deg(flf_alpha[0]), round)*3600, label="$\omega_{x\\alpha}^F$", alpha = 0.5);
fig4_ax1.grid(True)
fig4_ax1.legend(loc="best")
#
fig4_ax2.set_title("Оценка скоростей дрейфа гироскопа Y")
fig4_ax2.set_xlabel("мин")
fig4_ax2.set_ylabel("град/час")
fig4_ax2.plot(time, np.round(np.rad2deg(data["d_omega_y"]), round)*3600, label="$\delta\omega_y$");
# fig4_ax2.plot(time, np.round(np.rad2deg(mean[1]), round)*3600, label="$E(\omega_y^F)$");
# fig4_ax2.plot(time, np.round(np.rad2deg(flf[1]), round)*3600, label="$\omega_y^F$", alpha = 0.5, linestyle = "--", marker = "*");
fig4_ax2.plot(time, np.round(np.rad2deg(flf_alpha[1]), round)*3600, label="$\omega_{y\\alpha}^F$", alpha = 0.5);
fig4_ax2.grid(True)
fig4_ax2.legend(loc="best")

plt.show()
