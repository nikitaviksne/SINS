import numpy as np
import csv
import matplotlib.pyplot as plt

dim_x = 6;
dim_z = 2
x0 = np.zeros(dim_x)#начальное значение вектора состояния
# x0[:2] = Vins[:2]
g = 9.80665
U = np.deg2rad(15)/3600 #радианы в секунды
R_e = 6400*10**3 #(м) экваториальный радиус
a =	6378245 #большая полуось
b = 6356856
e = np.sqrt(1-b**2/a**2) #эксцентрисетет
height = 0;#высота
beta_ff = 6e-8
a_ff = 1e-4
freq = 100 #in Hz
dt = 1/freq # шаг 
t_model = (90-5)*60 #секунды
#N = int(t_nav/dt)
gyro_drift = np.deg2rad(0.05)/3600 #рад/с
acc_drift = 1e-3*g
Xplot = [[] for i in range(dim_x)] # для складирования оценненых значений вектора состояния
X = np.zeros((dim_x, 1)) # для складирования оценненых значений вектора состояния только на 1 такте
X[2, 0] = acc_drift/g; X[3, 0] = -acc_drift/g; # начальные ошибки ориентации (невыставка)

allow_bias_acc = True;
allow_rand_acc = True;
allow_rand_V = True;
allow_rand_Coo = True;
allow_bias_gyr = False
psi = np.deg2rad(90)

path = "/home/nikita_viksne/InertialNavigation/Data_files/data_acc_veloc_30_heading_90_freq_400_V_coo_gps.csv" #путь к файлу с моделированными показаниями чувствитлеьныхэлементов

Gyro = np.array([[gyro_drift], 
                 [gyro_drift]])
B_gyro = allow_bias_gyr*np.array([[0, 0], # Ve
                    [0, 0], # Vn
                    [0, 0], #Phi_e
                    [0, 0], #phi_n
                    [1, 0],#d_om_x
                    [0, 1]]) # d_om_y

with open(path, newline='') as file:
    iter = 0;
    reader = csv.reader(file, delimiter = ";")
    for _ in range(5*60*400): #пропускаем выставку
        next(reader) #костыль, но лучше чем file.seek(5*60*400), так как я не знаю точной длины каждой строки (количество символов разное)
    for row in reader:
        V_ideal = np.array(row[-9:-7], dtype = float) # иделаьные скорости от gps, но пусть будут называться 
        Vins = V_ideal + allow_rand_V*np.array(row[-7:-5], dtype = float) # иделаьные скорости от gps, но пусть будут называться ins
        Cooins = np.array(row[-5:-3], dtype = float) + allow_rand_Coo*np.array(row[-3:-1], dtype = float)
        A = np.array(row[0:3], dtype = float) + allow_bias_acc*np.array(row[6:9], dtype = float) + allow_rand_acc*np.array(row[12:15], dtype=float)
        omega_s = np.array([-Vins[0]/R_e, Vins[0]/(R_e*np.cos(Cooins[0])), Vins[0]/R_e*np.tan(Cooins[0])])

        A_m = np.array([[Vins[1]/(R_e)*np.tan(Cooins[0]), Vins[0]/(R_e)*np.tan(Cooins[0]) + 2*U*np.sin(Cooins[0]), 0, -A[2], 0, 0], # Delta dot V_ox
                    [-2*(Vins[0]/(R_e)*np.tan(Cooins[0]) + U*np.sin(Cooins[0])), 0, A[2], 0, 0, 0], # Delta dot V_oy
                    # [2*(V_ox/(R_e+height) + U*np.cos(Cooins[0])), 2*V_oy/(R_e+height), 0, -A[1], A[0], 0, 0], # Delta dot Vins[2]
                    [0, -1/(R_e), 0, omega_s[2], -np.cos(psi), -(np.sin(psi))], # Phi_ox
                    [1/(R_e), 0, - omega_s[2], 0, -(-np.sin(psi)), -(np.cos(psi))], # Phi_oy
                    # [np.tan(Cooins[0])/(R_e+height), 0, 0, omega_s[1], -omega_s[0], 0, 0, 0], # Delta dot Phi_oz
                    [0, 0, 0, 0, 0, 0], # Delta omega_x
                    [0, 0, 0, 0, 0, 0]]) # Delta omega_y
        X = A_m @ X + B_gyro@Gyro
        with open("./Data_files/Modelling_errors.csv", 'at') as outfile:
            outfile.write(f"{iter},")
            for i in range(3):
                outfile.write(f"{A[i]},")
            for i in range(3):
                outfile.write(f"{omega_s[i]},")
            for i in range(2):
                outfile.write(f"{V_ideal[i] + X[i,0]},") #скорости ИНС это идеал GPS + ошибка
            outfile.write(f"{psi},")
            outfile.write(f"0,") # pitch
            outfile.write(f"0,") # roll
            for i in range(2):
                outfile.write(f"{Cooins[i]},")
            for i in range(2):
                outfile.write(f"{Vins[i]},")
            for i in range(dim_x):
                outfile.write(f"{X[i,0]},")
            outfile.write("\n")
        
        for i in range(dim_x):
            Xplot[i].append(X[i])
        
        iter += 1 

        for _ in range(3):
            try:
                next(reader)    #пропускаю 3 такта, так как в файле частота 400 Гц, а я моделирую с частотой 100 Гц
            except StopIteration:
                break;


time = np.linspace(0, len(Xplot[0])/freq/60, len(Xplot[0]))

fig1, (fig1_ax1, fig1_ax2) = plt.subplots(2, 1)
fig1.suptitle("Скорости")
fig1_ax1.set_title("Восточная скорость");
fig1_ax1.set_xlabel("мин")
fig1_ax1.set_ylabel("м/с")
fig1_ax1.plot(time, Xplot[0])
fig1_ax1.grid(True)
#
fig1_ax2.set_title("Северная скорость");
fig1_ax2.set_xlabel("мин")
fig1_ax2.set_ylabel("м/с")
fig1_ax2.plot(time, Xplot[1])
fig1_ax2.grid(True)

'''
fig2, (fig2_ax1) = plt.subplots(1, 1)
fig2_ax1.plot(time, np.rad2deg(Xplot[1])*60)
'''

fig3, (fig3_ax1) = plt.subplots(1, 1)
fig3.suptitle("Скорость дрейфа гироскопа")
fig3_ax1.set_xlabel("мин")
fig3_ax1.set_ylabel("град/ч")
fig3_ax1.plot(time, np.rad2deg(Xplot[-1])*3600)
fig3_ax1.axhline(y = 0.05, xmin=time[0], color = "#FF0000", label = "Модель")
fig3_ax1.grid(True)

'''
fig4, (fig4_ax1) = plt.subplots(1, 1)
fig4_ax1.plot(time, np.rad2deg(Orient[2]))
'''

plt.show()
