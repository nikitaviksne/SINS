'''
Моделирует показания чувствительных элементов при движении с заданными кглами ориентации (углами Эйлера-Крылова)
и с заданными скоростями. Все сохраняет в csv файл, разделенный пробелами (так C++ лучше читает)
'''
import numpy as np
import struct

def matrix_o_b(H, R, P):
    '''Функция формирует матрицу перехода (МНК) из опорной систему координат (географической) в связанную'''
    '''
    H - Heading (курс)
    R - Roll (крен)
    P - Pitch (тангаж)
    '''
    C = np.array([
        [np.cos(R) * np.cos(H) + np.sin(R)*np.sin(H)* np.sin(P), np.sin(R)*np.cos(H)*np.sin(P) - np.cos(R)*np.sin(H), -np.sin(R) * np.cos(P)],
        [np.cos(P)*np.sin(H), np.cos(H)*np.cos(P), np.sin(P)],
        [np.sin(R)*np.cos(H) - np.cos(R)*np.sin(H)*np.sin(P), -np.sin(R)*np.sin(H) - np.cos(R)*np.cos(H)*np.sin(P), np.cos(R)*np.cos(P)]
    ])
    return C

def make_portable_speed(phi):
#     portable -- переносная. Формирование переносной угловой скорости (т.е. скорости Земли)
#     В осях опорной системы координат
    Om_p = np.array([0, U*np.cos(phi), U*np.sin(phi)]);
    return Om_p


def cross_sim(A):
    '''Кососиметрическая матрица для угловых скоростей'''
    return np.array([0,-A[2], A[1],A[2],0,-A[0], -A[1],A[0],0]).reshape(3,3)

U = np.deg2rad(15)/3600;
g = 9.81;
Re = 6400e3; # радиус Земли в метрах
a = 6378245;
b = 6356856;
e=np.sqrt(1 - b*b/a/a);

freq = 400; # частота (измерений за 1 секунду)
t_nav = 90*60 # в секундах


'''
C_b_n = np.array([[],
                  [],
                  []])
'''

Vabs = 0 # модуль конечной линейной скорости, м/с
phi0 = np.deg2rad(0)

'''задаем ориентацию объекта'''
heading = np.deg2rad(0)
roll = np.deg2rad(0);
pitch = np.deg2rad(0);

file_name = f"data_acc_veloc_{Vabs}_heading_{int(np.rad2deg(heading))}_eqautor_grinvich_freq_400_gps.csv"
C_n_b = matrix_o_b(heading, roll, pitch)

'''Систематические дрейфы'''
bias_acc = 1e-4;
bias_gyr = np.deg2rad(0.05)/3600;
'''Случайные дрейфы'''
T_k_a = 1
beta_acc = 1/T_k_a
std_acc = 0.1*g*1e-3
T_k_g = 2
beta_gyr = 1/T_k_g
std_gyr = np.deg2rad(0.02)/3600 #rad in sec

num_samples = t_nav*freq;
mean = 0;

white_noise_acc= np.random.normal(mean, 1, size=(3, num_samples))

white_noise_gyr = np.random.normal(mean, 1, size=(3, num_samples))

colour_noise_acc = np.zeros((3, num_samples))
colour_noise_gyr = np.zeros((3, num_samples))


sko_gnss_n = 0.05 #in meter in seconds
sko_gnss_d = sko_gnss_n*np.sqrt(freq)
# sko_gnss_d = 0.2
sko_gnss_pos = 0.2/Re*np.sqrt(freq)

'''формирующий фильтр для экспоненциально косинусоидальной корреляционной функции'''
for i in range(1,num_samples):
    for jj in range(3):
        #акселерометры
        colour_noise_acc[jj][i] = colour_noise_acc[jj][i-1] * (1 - beta_acc/freq) + std_acc*np.sqrt(2*beta_acc/freq)*white_noise_acc[jj][i-1]
        #гироскопы
        colour_noise_gyr[jj][i] = colour_noise_gyr[jj][i-1]*(1-beta_gyr/freq) + std_gyr*np.sqrt(2*beta_gyr/freq)*white_noise_gyr[jj][i-1]


'''Моделируем измерения в связанной системе координат и записываем в файл'''
time_to_alignment = 5*60; # время выставки в секундах
Ve = np.linspace(Vabs*np.sin(heading), Vabs*np.sin(heading), (t_nav - time_to_alignment)*freq)
Vn = np.linspace(Vabs*np.cos(heading), Vabs*np.cos(heading), (t_nav - time_to_alignment)*freq)

Ve = np.concatenate((np.array([0 for _ in range(time_to_alignment*freq)]), Ve))
Vn = np.concatenate((np.array([0 for _ in range(time_to_alignment*freq)]), Vn))

# Изменение высоты полета
Height = np.linspace(0, 0, (t_nav - time_to_alignment)*freq)
Height = np.concatenate((np.zeros(time_to_alignment*freq), Height))
# Изменение широты
phi = [phi0]

for iii in range(1, num_samples): # цикл от 1 (практически с начала)
    phi.append(phi[-1] + Vn[iii]/(Re + Height[iii])/freq)


round = 12; # количество знаков после запятой, с какой округлять и выводить в файл
with (open(f"./Data_files/{file_name}", "wt") as file):
    '''
    последовательность данных
    Abx Aby Abz Ombx Omby Omz biasAbx biasAby biasAbz biasOmbx biasOmby biasOmbz randAbx randAby randAbz randOmbx randOmby randOmbz Vgps_x Vgps_y randVgps_x randVgps_y
    '''
    s = struct.Struct("<22d");
    for itr in range(num_samples - 1):
        Rlambda = Re# / np.sqrt(1 - e **2 * np.sin(phi[itr]) **2);
        Rphi = Re# * (1 - e**2) / (np.sqrt(1 - e**2 * np.sin(phi[itr]) **2) * (1 - e**2 * np.sin(phi[itr]) **2));

        A_o = np.array([0.0, 0.0, g])

        Om_e = make_portable_speed(phi[itr]); # формирование скорости вращения Земли на широте

        # Дифференцирование скоростей и прибавление к ускорению
        if (itr>time_to_alignment*freq):
            dA = np.array([Ve[itr + 1] - Ve[itr], Vn[itr + 1] - Vn[itr],
                           Height[itr +1] - Height[itr]])
            # Изменение угловых скоростей от движения
            dOm_or = np.array([-Vn[itr] / (Rphi + Height[itr]), #
                              Ve[itr] / ( (Rlambda + Height[itr]) ), #* np.cos(phi[itr] )
                              Ve[itr] * np.tan(phi[itr]) / ((Rlambda + Height[itr]) ) ]) # r -- realtive (относительная)
            '''
            Ускорение Кориолиса равно удвоенному векторному произведению угловой
            скорости переносного движения на относительную скорость точки (видимо имеется в виду линейную скорость)
            При движении по поверхности Земли, переностная скорость --- скорость вращения земли (0, U*cos, U*sin), относительная -- скорость относительно наблюдателя на Земле.
            '''
            Coriolise = np.cross(2*(Om_e), np.array([Ve[itr], Vn[itr], 0])); #+  np.cross(Om_e, np.array([Ve[itr], Vn[itr], 0]))
        else:
            dA = np.zeros(3)
            dOm_or = np.zeros(3) # r -- realtive (относительная)
            Coriolise = np.zeros(3);

        Om_b = C_n_b @ (Om_e + dOm_or);
        A_b = C_n_b @ (A_o + dA + Coriolise);

        '''Запись в файл в текстовом виде'''
        # '''
        for i in range(3):
            file.write(f"{np.round(A_b[i], round)};")
        for i in range(3):
            file.write(f"{np.round(Om_b[i], round)};")

        for i in range(3):
            file.write(f"{np.round(bias_acc, round)};") # Смещение нулей акселерометров
        for i in range(3):
            file.write(f"{np.round(bias_gyr, round)};") # Смещение нулей гироскопов в рад/с
        for i in range(3):
            file.write(f"{np.round(colour_noise_acc[i][itr], round)};") # Случайные погрешности акселерометров
        for i in range(3):
            file.write(f"{np.round(colour_noise_gyr[i][itr], round)};") # Случайные погрешности гироскопов

        file.write(f"{np.round(Ve[itr], round)};") # Восточная скорость (идеал, будет воспринята как GPS)

        file.write(f"{np.round(Vn[itr], round)};") # Северная скорость (идеал, будет воспринята как GPS)

        file.write(f"{np.round(np.random.normal(loc = 0, scale = sko_gnss_d), round)};") # Один белый гауссовский шум для скорости GPS
        file.write(f"{np.round(np.random.normal(loc = 0, scale = sko_gnss_d), round)};") # Второй белый гауссовский шум для скорости GPS

        file.write("\n")
        # '''

        '''
        data = s.pack(*A_b, *Om_b, bias_acc, bias_acc, bias_acc, bias_gyr, bias_gyr, bias_gyr, *colour_noise_acc[:, itr], *colour_noise_gyr[:, itr], Ve[itr], Vn[itr], np.random.normal(loc = 0, scale = sko_gnss_d), np.random.normal(loc = 0, scale = sko_gnss_d) ) # формирование пакета (структуры)

        file.write(data); # запись этой структуры в файл
        # '''

print("Generation done")
