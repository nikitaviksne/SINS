import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys #для приема имени файла для построения через аргумент команды

path_nav_sol = sys.argv[1]
path_est = sys.argv[2]
# data =pd.read_csv("~/InertialNavigation/data/"+relative_path, delimiter=";");
data = pd.read_csv(path_nav_sol, delimiter=";");
estimations = pd.read_csv(path_est, delimiter=",");


time = np.linspace(0, (data.iloc[:, 0].size - 1)/100/60, data.iloc[:, 0].size) #в минутах /100/60

round = 100

'''Ошибки линейный скоростей'''
fig1, (fig1_ax1, fig1_ax2) = plt.subplots(1, 2)
# Восточная составляющая
fig1.suptitle("Ошибки скоростей")
fig1_ax1.set_title("Ошибка восточной скорости");
fig1_ax1.set_xlabel("мин")
fig1_ax1.set_ylabel("м/с")
fig1_ax1.plot(time, np.round(data["errVe"], round), label="$\Delta$ Ve");
fig1_ax1.plot(time, np.round(estimations["Ve"], round), linestyle='--', label="$\hat{\Delta Ve}$");
fig1_ax1.legend(loc="best")
fig1_ax1.grid(True)
# Северная составляющая
fig1_ax2.set_title("Ошибка северной скорости");
fig1_ax2.set_xlabel("мин")
fig1_ax2.set_ylabel("м/с")
fig1_ax2.plot(time, np.round(data["errVn"],round), label="$\Delta$ Vn");
fig1_ax2.plot(time, np.round(estimations["Vn"],round), linestyle='--', label="$\hat{\Delta Vn}$");
fig1_ax2.legend(loc="best")
fig1_ax2.grid(True)

'''Углы ориентации'''
fig2, (fig2_ax1, fig2_ax2, fig2_ax3) = plt.subplots(3, 1)
fig2.suptitle("Углы")
# Угол курса
fig2_ax1.set_title("Угол курса");
fig2_ax1.set_xlabel("мин")
fig2_ax1.set_ylabel("Град.")
fig2_ax1.plot(time, np.round(np.rad2deg(data["Heading"]),round), label="Курс");
fig2_ax1.plot(time, np.round(np.rad2deg(estimations["d_Psi"]),round), label="$\hat{\Delta\Psi}$");
fig2_ax1.legend(loc="best")
fig2_ax1.grid(True)
# Угол крена
fig2_ax2.set_title("Угол крена");
fig2_ax2.set_xlabel("Мин")
fig2_ax2.set_ylabel("Угл. мин")
fig2_ax2.plot(time, np.round(np.rad2deg(data["Roll"])*60,round), label="$\gamma$");
fig2_ax2.plot(time, np.round(np.rad2deg(estimations["d_Roll"])*60,round), label="$\hat{\Delta\gamma}$");
fig2_ax2.legend(loc="best")
fig2_ax2.grid(True)
# Угол тангажа
fig2_ax3.set_title("Угол тангажа");
fig2_ax3.set_xlabel("Мин")
fig2_ax3.set_ylabel("Угл. мин")
fig2_ax3.plot(time, np.round(np.rad2deg(data["Pitch"])*60,round), label="$\\theta$");
fig2_ax3.plot(time, np.round(np.rad2deg(estimations["d_Pitch"])*60,round), label="$\hat{\Delta\\theta}$");
fig2_ax3.legend(loc="best")
fig2_ax3.grid(True)

'''Осредняем скорости дрейфов гироскопов и сглаживаем ФНЧ с T=100'''

Tf = 1e6;
alpha_f = 1e-5;
mean = [[0],[0]] # среднее значение ФНЧ
flf = [[0], [0]] # filter lower frequency
flf_alpha = [[0], [0]] # filter lower frequency с альфой
for i, dx in enumerate(estimations["d_omega_x"]):
    #MeanAb[i] = (Ldoub) MeanAb[i] + (Ab[i] - MeanAb[i]) / (iter + 1); //(iter * MeanAb[i] + Ab[i])/(iter + 1);
    flf[0].append((dx + Tf/100*flf[0][-1])/(1 + Tf/100))
    flf[1].append((estimations["d_omega_y"][i+1] + Tf/100*flf[1][-1])/(1 + Tf/100))
    # Фильтр с альфой
    flf_alpha[0].append(flf_alpha[0][-1]*(1-alpha_f) + dx*alpha_f);
    flf_alpha[1].append(flf_alpha[1][-1]*(1-alpha_f) + estimations["d_omega_y"][i+1]*alpha_f);
    mean[0].append(mean[0][i] + (flf[0][i+1] - mean[0][i])/(i + 1));
    mean[1].append(mean[1][i] + (flf[1][i+1] - mean[1][i])/(i + 1));
del(mean[0][0])
del(mean[1][0])
del(flf[0][0])
del(flf[1][0])
del(flf_alpha[0][0])
del(flf_alpha[1][0])

'''Оценки дрейфов гироскопов'''
fig3, (fig3_ax1, fig3_ax2) = plt.subplots(2, 1)
# Восточное направление
fig3.suptitle("Оценка скоростей дрейфов гироскопов")
fig3_ax1.set_title("Оценка скоростей дрейфа гироскопа канала X")
fig3_ax1.set_xlabel("Мин")
fig3_ax1.set_ylabel("Град/час")
'''
fig3_ax1.axhline(y = 0.05, xmin=0, color = "#FF0000", label = "Модель, body")
fig3_ax1.axhline(y = 0.05*np.cos(np.deg2rad(50)) + 0.05*np.sin(np.deg2rad(50)), xmin=0, color = "#db0db9", linestyle = "-.", label = "Модель, psi = 50$\circ$")
fig3_ax1.axhline(y = -0.05*np.sin(np.deg2rad(50)) + 0.05*np.cos(np.deg2rad(50)), xmin=0, color = "#db0db9",  linestyle = "-.", label = "Модель, psi = 50$\circ$")
fig3_ax1.axhline(y = 0.05*np.cos(np.deg2rad(50+90)) + 0.05*np.sin(np.deg2rad(50+90)), xmin=0, color = "#0ddb60", linestyle = "--", label = "Модель, psi = 140$\circ$")
fig3_ax1.axhline(y = -0.05*np.sin(np.deg2rad(50+90)) + 0.05*np.cos(np.deg2rad(50+90)), xmin=0, color = "#0ddb60",  linestyle = "--", label = "Модель, psi = 140$\circ$")
'''
fig3_ax1.plot(time, np.round(np.rad2deg(estimations["d_omega_x"]), round)*3600, label="$\delta\omega_x$");
# fig3_ax1.plot(time, np.round(np.rad2deg(mean[0]), round)*3600, label="$E(\omega_x^F)$");
# fig3_ax1.plot(time, np.round(np.rad2deg(flf[0]), round)*3600, label="$\omega_x^F$", alpha = 0.5, linestyle = "--", marker = "*");
fig3_ax1.plot(time, np.round(np.rad2deg(flf_alpha[0]), round)*3600, label="$\omega_{x}^{F\\alpha}$");
fig3_ax1.grid(True)
fig3_ax1.legend(loc="best")
#
fig3_ax2.set_title("Оценка скоростей дрейфа гироскопа канала Y")
fig3_ax2.set_xlabel("Мин")
fig3_ax2.set_ylabel("Град/час")
'''
fig3_ax2.axhline(y = 0.05, xmin=0, color = "#FF0000", label = "Модель, body")
fig3_ax2.axhline(y = 0.05*np.cos(np.deg2rad(50)) + 0.05*np.sin(np.deg2rad(50)), xmin=0, color = "#db0db9", linestyle = "-.", label = "Модель, psi = 50$\circ$")
fig3_ax2.axhline(y = -0.05*np.sin(np.deg2rad(50)) + 0.05*np.cos(np.deg2rad(50)), xmin=0, color = "#db0db9",  linestyle = "-.", label = "Модель, psi = 50$\circ$")
fig3_ax2.axhline(y = 0.05*np.cos(np.deg2rad(50+90)) + 0.05*np.sin(np.deg2rad(50+90)), xmin=0, color = "#0ddb60", linestyle = "--", label = "Модель, psi = 140$\circ$")
fig3_ax2.axhline(y = -0.05*np.sin(np.deg2rad(50+90)) + 0.05*np.cos(np.deg2rad(50+90)), xmin=0, color = "#0ddb60",  linestyle = "--", label = "Модель, psi = 140$\circ$")
'''
fig3_ax2.plot(time, np.round(np.rad2deg(estimations["d_omega_y"]), round)*3600, label="$\delta\omega_y$");
# fig3_ax2.plot(time, np.round(np.rad2deg(mean[1]), round)*3600, label="$E(\omega_y^F)$");
# fig3_ax2.plot(time, np.round(np.rad2deg(flf[1]), round)*3600, label="$\omega_y^F$", alpha = 0.5, linestyle = "--", marker = "*");
fig3_ax2.plot(time, np.round(np.rad2deg(flf_alpha[1]), round)*3600, label="$\omega_{y}^{F\\alpha}$");
fig3_ax2.grid(True)
fig3_ax2.legend(loc="best")


plt.show()

exit();

'''Оценки ориентации'''
fig10, (fig10_ax1, fig10_ax2) = plt.subplots(2, 1)
# Угол Phi_x
fig10_ax1.set_title("Оценка угла $\Phi_x$");
fig10_ax1.set_xlabel("Мин")
fig10_ax1.set_ylabel("Угл. мин")
fig10_ax1.plot(time, np.round(np.rad2deg(estimations["Phi_e"])*60,round), label="$\hat{\Phi_e}$");
fig10_ax1.legend(loc="best")
fig10_ax1.grid(True)
# Угол phi_y
fig10_ax2.set_title("Оценка угла $\Phi_n$");
fig10_ax2.set_xlabel("Мин")
fig10_ax2.set_ylabel("Угл. мин")
fig10_ax2.plot(time, np.round(np.rad2deg(estimations["Phi_n"])*60,round), label="$\hat{\Phi_n}$");
fig10_ax2.legend(loc="best")
fig10_ax2.grid(True)



'''Ошибки по координатам'''
fig4, (fig4_ax1, fig4_ax2) = plt.subplots(1, 2)
# Восточное направление
fig4.suptitle("Ошибка по координатам")
fig4_ax1.set_title("Восточное направление");
fig4_ax1.set_xlabel("Мин")
fig4_ax1.set_ylabel("М")
fig4_ax1.plot(time, np.round(data["d_E"], round), label="dE");
fig4_ax1.legend(loc="best")
fig4_ax1.grid(True)
# Северное направление
fig4_ax2.set_title("Северное направление");
fig4_ax2.set_xlabel("Мин")
fig4_ax2.set_ylabel("М")
fig4_ax2.plot(time, np.round(data["d_N"],round), label="dN");
fig4_ax2.legend(loc="best")
fig4_ax2.grid(True)

'''Траектория'''
fig5, (fig5_ax1) = plt.subplots(1, 1)
# Восточное направление
fig5.suptitle("Траектория")
fig5_ax1.set_xlabel("град")
fig5_ax1.set_ylabel("град")
fig5_ax1.plot(np.round(np.rad2deg(data["Lambda"]), round), np.round(np.rad2deg(data["Phi"]), round) );
fig5_ax1.grid(True)

plt.show()
