import csv
import numpy as np
from RAdaptiveKalman import *
import matplotlib.pyplot as plt


allow_bias_acc = [True]*3
allow_bias_gyr = [True]*3
allow_rand_acc = [True]*3
allow_rand_gyr = [True]*3
allow_rand_Vgps = [True]*2
allow_rand_Coogps = [False]*2

path = "/home/nikita_viksne/Документы/Python/DevelopKalman/DataFiles/For_Kalman_freq_400.csv"
#чтение файла
with open(path, newline='') as file:
    iter = 0;
    reader = csv.reader(file, delimiter = ",")
    print(*enumerate(next(reader))) # выбрасываем строку заголовков
    
    row = next(reader)
    Vgps = [0]*3 # скорости с gps
    Vins = [0]*3 # скорости с INS (Z пока нулевая)
    Cooins = [0]*3 # координаты с ins
    for i in range(2):
        Vins[i] = float(row[7 + i])
    
    dim_x = 6;
    dim_z = 2
    x0 = np.zeros(dim_x)#начальное значение вектора состояния
    # x0[:2] = Vins[:2]
    filter = R_Adaptive_Kalman(dim_x, dim_z, x=x0)
    X = [[] for i in range(dim_x)] # для складирования оценненых значений вектора состояния
    g = 9.80665
    U = np.deg2rad(15)/3600 #радианы в секунды
    R_e = 6400*10**3 #(м) экваториальный радиус
    a =	6378245 #большая полуось
    b = 6356856
    e = np.sqrt(1-b**2/a**2) #эксцентрисетет
    height = 0;#высота
    beta_ff = 6e-8
    a_ff = 1e-4
    freq = 100
    dt = 1/freq # шаг 
    #t_nav = 90*60 #секунды
    #N = int(t_nav/dt)
    Orient = [[] for i in range(3)];
    '''
    Параметры модели
    '''

    
    # G = np.array([[0], [0], [a_ff*np.sqrt(2*beta_ff)]])
    # G = np.zeros((3,3)); G[2,:] = np.array([0, 0, a_ff*np.sqrt(2*beta_ff)*dt])
    G = np.array([[0, 0, 0, 0, 1*dt, 0], 
                [0, 0, 0, 0, 0, 1*dt]]).T
    H = np.array([[1, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0]])
    Q = np.eye(2)*1e-17
    
    filter.G = G
    filter.H = H
    filter.Q = Q

    for row in reader: #в row список строк 
        
        A = [0]*3 # ускорения с акселерометров
        w = [0]*3 #угловые скорости с гироскопов
        Vgps = [0]*3 # скорости с gps
        Vins = [0]*3 # скорости с INS (Z пока нулевая)
        Cooins = [0]*3 # координаты с ins
        Coogps = [0]*3 # координаты с gps
        Cbn = np.zeros((3,3)) #матрица перехода из свзяанной в опорную с.к
         
        omega_s = [0]*3 # скорости с INS (Z пока нулевая)
        #Timestamp,a_ll_x,a_ll_y,a_ll_z,omega_s_x,omega_s_y,omega_s_z,Ve_ins,Vn_ins,heading,roll,pitch,latitude_ins,longitude_ins,Ve_gps,Vn_gps,C00,C01,C02,C10,C11,C12,C20,C21,C22,
        heading = float(row[9])
        Orient[2].append(heading)
        for i in range(3): 
            A[i] += float(row[1 + i]) # формирование измерений с инерциальных ч.э.
            #w[i] += float(row[4 + i]) # формирование измерений с инерциальных ч.э.
            omega_s[i] = float(row[4 + i]) # абсолютные угловые скорости в опорной с.к
        for i in range(2):
            Vins[i] = float(row[7 + i])
            Cooins[i] = float(row[12 + i])
            Vgps[i] += float(row[14 + i]) # формирование измерений с неинерциальных ч.э. (ГНСС)
        for iii in range(3): # по строкам
            for jjj in range(3): # по столбцам
                Cbn[iii,jjj] = float(row[16 + 3*iii+jjj])

        '''    
        for i in range(2): # формирование измерений с неинерциальных ч.э. (ГНСС)
            Vgps[i] += float(row[7 + i])
            #Coogps[i] += 0*float(row[22 + i]) + allow_rand_Coogps[i] * float(row[24 + i])
        '''
        A_m = np.array([[Vins[1]/(R_e)*np.tan(Cooins[0]), Vins[0]/(R_e)*np.tan(Cooins[0]) + 2*U*np.sin(Cooins[0]), 0, -A[2], 0, 0], # Delta dot V_ox
                    [-2*(Vins[0]/(R_e)*np.tan(Cooins[0]) + U*np.sin(Cooins[0])), 0, A[2], 0, 0, 0], # Delta dot V_oy
                    # [2*(V_ox/(R_e+height) + U*np.cos(Cooins[0])), 2*V_oy/(R_e+height), 0, -A[1], A[0], 0, 0], # Delta dot Vins[2]
                    [0, -1/(R_e), 0, omega_s[2], -Cbn[0,0], -Cbn[0,1]], # Phi_ox
                    [1/(R_e), 0, - omega_s[2], 0, -Cbn[1,0], -Cbn[1,1]], # Phi_oy
                    # [np.tan(Cooins[0])/(R_e+height), 0, 0, omega_s[1], -omega_s[0], 0, 0, 0], # Delta dot Phi_oz
                    [0, 0, 0, 0, 0, 0], # Delta omega_x
                    [0, 0, 0, 0, 0, 0]]) # Delta omega_y
                  #   [0, 0, 0, 0, 0, 0, 0, 0, -beta]]) # Delta omega_z
        #A_m = np.array([[0, -A[2], 0], [1/R_e, 0, 1], [0, 0, -beta_ff]])
        Phi = np.eye(dim_x) + A_m*dt
        filter.Phi = Phi
        
        z = np.array([Vins[i] - Vgps[i] for i in range(dim_z)]) # 
        filter.predict()
        filter.update(z)

        #Debug
        '''
        print(f"iter = {iter}")
        print(f"x_k_1 = {filter.x_k_1}")
        print(f"P_apriori = {filter.P_apriori}")
        print(f"v_k = {filter.v_k}")
        print(f"C_k = {filter.C_k}")
        print(f"R_k = {filter.R_k}")
        print(f"k = {filter.k}")
        print(f"P_aposteriori = {filter.P_aposteriori}")
        print(f"x = {filter.x}")
        print(f"k_iter = {filter.k_iter}")
        # print(f"P_apriori = {filter.P_apriori}")
        # '''
        
        for i in range(dim_x):
            X[i].append(filter.x[i])
        
        '''
        iter += 1
        if iter > 2: #20*60*freq:
            break;
        # '''

        '''
        print(f"A = {A}")
        print(f"w = {np.rad2deg(w)*3600}")
        print(f"Vgps = {Vgps}")
        print(f"Coogps = {Coogps}")
        iter += 1
        if iter>10:
            break;
        '''

time = np.linspace(0, len(X[0])/freq/60, len(X[0]))

fig1, (fig1_ax1, fig1_ax2) = plt.subplots(2, 1)
fig1.suptitle("Скорости")
fig1_ax1.set_title("Восточная скорость");
fig1_ax1.set_xlabel("мин")
fig1_ax1.set_ylabel("м/с")
fig1_ax1.plot(time, X[0])
fig1_ax1.grid(True)
#
fig1_ax2.set_title("Северная скорость");
fig1_ax2.set_xlabel("мин")
fig1_ax2.set_ylabel("м/с")
fig1_ax2.plot(time, X[1])
fig1_ax2.grid(True)

'''
fig2, (fig2_ax1) = plt.subplots(1, 1)
fig2_ax1.plot(time, np.rad2deg(X[1])*60)
'''

fig3, (fig3_ax1) = plt.subplots(1, 1)
fig3.suptitle("Скорость дрейфа гироскопа")
fig3_ax1.set_xlabel("мин")
fig3_ax1.set_ylabel("град/ч")
fig3_ax1.plot(time, np.rad2deg(X[-1])*3600)
fig3_ax1.axhline(y = 0.05, xmin=time[0], color = "#FF0000", label = "Модель")
fig3_ax1.grid(True)

'''
fig4, (fig4_ax1) = plt.subplots(1, 1)
fig4_ax1.plot(time, np.rad2deg(Orient[2]))
'''

plt.show()
