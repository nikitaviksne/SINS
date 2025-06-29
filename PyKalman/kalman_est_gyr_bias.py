import numpy as np
import pickle
import pandas as pd
from Kalman import Adaptive_Kalman, Kalman
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt


'''
Попробуем сделать макет для фильтра калмана для оценки дрейфов по модели которую составил я
'''


U = np.radians(15)/3600
R = 6400*1e3
R_e = 6400*10**3 #(м) экваториальный радиус
g = 9.80665
# beta = 1*6e-10
beta = 2.5e-1
# beta = 2.5
# a_ff = 1e-4
tau = 1e-2
t_nav = 2*40*60
time_start_nav = int(5*60/tau) #время выставки в тактах
# N = int(t_nav/tau) -time_start_nav
# N = int(t_nav/tau) -time_start_nav



'''Открываю файл для чтения данных решения автономного алгорима'''
file_path_ins = "/home/nikita_viksne/InertialNavigation/data/For_Kalman.csv"#"./Data_files/For_Kalman.csv"
data = pd.read_csv(file_path_ins, delimiter=",", index_col = False);
N = len(data["Timestamp"]) # клоличество строк, задается вручную

'''создание файла для записи оценок дрейфов'''
file_estimation_path = "./Data_files/Saratov_estimation.csv"

allow_gyr_bias = 1 # включение смещение нуля ДУС

dim_x = 7
dim_z = 3
x_est = np.zeros((dim_x, N ))
omega_est = np.zeros((2, N ))

filter = Kalman(dim_x = dim_x, dim_z = dim_z, x = np.zeros(dim_x))
# filter = Adaptive_Kalman(dim_x = dim_x, dim_z = dim_z, x = np.zeros(dim_x))
filter.P_aposteriori = np.zeros((dim_x, dim_x))
filter.Q = np.eye(dim_z)*1e-16#*file_alg["h"]#; filter.Q[0,0] = 0
# filter.Q = np.array([1e-13])
# filter.G = np.array([0, 0, 0, 0, 0, 0, 1*tau, 1*tau, 1*tau]).reshape(dim_x,1)
filter.G = np.zeros((dim_x, dim_z)); filter.G[-2:,-2:] = np.eye(2)*1e-2 # 1e-2 -- период дискретизации (100 Гц)
# filter.H = np.array([1, 0, 0, 0, 0, 0, 0, 0]).reshape(1,8)
filter.H = np.zeros((dim_z,dim_x)); filter.H[0:dim_z,0:dim_z] = np.eye(dim_z)
I = np.eye(dim_x)
filter.R = np.eye(dim_z)/np.sqrt(0.01)*(0.05**2) #0.2 = scale in np.random.normal
time = []

# k_1 = file_alg["k_1"]
# k_2 = file_alg["k_2"]
nu = g/R
# for i in range(0, int(t_nav/tau)):
with open(file_estimation_path, "w") as est_file:# открываю файл для записи оценко дрейфов
    est_file.write("Timestamp,bias_est_bx,bias_est_by\n")#записываю заголовки
    for i in range(1, N-1):
      '''
      V_ox = gps_data['VeJPS'][i]
      V_oy = gps_data['VnJPS'][i]
      V_oz = gps_data['VyJPS'][i]
      # '''
    
      time.append(data["Timestamp"][i]) # в минутах
      """В качестве измерений попробую подавать ошибку по скорости между ИНС и GPS"""
      V_o_gps = np.array([data["Ve_gps"][i], data["Vn_gps"][i], 0])
      # V_o_gps = file_model["Velocity"][:, i] + np.random.normal(loc = 0, scale = 0.2, size = dim_z) # измерения с GPS истинная скорость с GPS
      V_ox = data["Ve_ins"][i] # выход с ИНС
      V_oy = data["Vn_ins"][i] # выход с ИНС
      V_ins = np.array([V_ox, V_oy, 0]) # массив numpy скорстей с инерциалки
      V_oz = 0
      height = 0
      # V_oz = 0
      # phi = file_alg['coord_auto'][0, i- file_model["time_to_align"]]
      # a_ox = file_alg['a_ll_auto'][0, i - file_model["time_to_align"]]
      # a_oy = file_alg['a_ll_auto'][1, i - file_model["time_to_align"]]
      # a_oz = file_alg['a_ll_auto'][2, i - file_model["time_to_align"]]
      # omega_ox = file_alg['omega_x'][i]
      # omega_oy = file_alg['omega_y'][i]
      # omega_oz = file_alg['omega_z'][i]
      # heading = file_model['heading_mas'][i]
      phi = data["latitude"][i]
      a_ox = data["a_ll_x"][i]
      a_oy = data["a_ll_y"][i]
      a_oz = data["a_ll_z"][i]
      omega_ox = data["omega_s_x"][i]
      omega_oy = data["omega_s_y"][i]
      omega_oz = data["omega_s_z"][i]
      heading = data["heading"][i]
      '''
      A = np.array([[V_oy/R*np.tan(phi) - V_oz/R, V_ox/R*np.tan(phi) + 2*U*np.sin(phi),  -(V_ox/R + 2*U*np.cos(phi)), 0,-a_oz, a_oy, 0, 0, 0],
                          [-2*V_ox/R*np.tan(phi) - 2*U*np.sin(phi), -V_oz/R, -V_oy/R, a_oz, 0, -a_ox, 0, 0, 0],
                          [2 * (V_ox/R + U * np.cos(phi)), 2 * V_oy/R, 0, -a_oy, a_ox, 0, 0, 0, 0],
                          [0, -1/R, 0, 0, omega_oz, -omega_oy, -np.cos(heading), -np.sin(heading), 0],
                          [1/R, 0, 0, -omega_oz, 0, omega_ox, np.sin(heading), -np.cos(heading), 0],
                          [np.tan(phi)/R, 0, 0, omega_oy, -omega_ox, 0, 0, 0, -1],
                          [0, 0, 0, 0, 0, 0, -beta, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, -beta, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, -beta]])
      # '''
      '''
      # модель из уравнений ошибок, упрощенная, с без курсового дрейфа и вертикального канала, размерность 6 
      A = np.array([[V_oy/(R+height)*np.tan(phi) - V_oz/(R+height), V_ox/(R+height)*np.tan(phi) + 2*U*np.sin(phi), 0, -a_oz, 0, 0], # Delta dot V_ox
                    [-2*(V_ox/(R+height)*np.tan(phi) + U*np.sin(phi)), -V_oz/(R+height), a_oz, 0, 0, 0], # Delta dot V_oy
                    # [2*(V_ox/(R+height) + U*np.cos(phi)), 2*V_oy/(R+height), 0, -a_oy, a_ox, 0, 0], # Delta dot V_oz
                    [0, -1/(R+height), 0, omega_oz, -np.cos(heading), -np.sin(heading)], # Phi_ox
                    [1/(R+height), 0, - omega_oz, 0, np.sin(heading), -np.cos(heading)], # Phi_oy
                  #   [np.tan(phi)/(R+height), 0, 0, omega_oy, -omega_ox, 0, 0, 0], # Delta dot Phi_oz
                    [0, 0, 0, 0, 0, 0], # Delta omega_x
                    [0, 0, 0, 0, 0, 0]]) # Delta omega_y
                  #   [0, 0, 0, 0, 0, 0, 0, 0, -beta]]) # Delta omega_z
      # '''
      # '''
      # модель из уравнений ошибок, упрощенная, с без курсового дрейфа, размерность 7
      A = np.array([[V_oy/(R+height)*np.tan(phi) - V_oz/(R+height), V_ox/(R+height)*np.tan(phi) + 2*U*np.sin(phi), -(V_ox/(R+height) + 2*U*np.cos(phi)), 0, -a_oz, 0, 0], # Delta dot V_ox
                    [-2*(V_ox/(R+height)*np.tan(phi) + U*np.sin(phi)), -V_oz/(R+height), -V_oy/(R+height), a_oz, 0, 0, 0], # Delta dot V_oy
                    [2*(V_ox/(R+height) + U*np.cos(phi)), 2*V_oy/(R+height), 0, -a_oy, a_ox, 0, 0], # Delta dot V_oz
                    [0, -1/(R+height), 0, 0, omega_oz, -np.cos(heading), -np.sin(heading)], # Phi_ox
                    [1/(R+height), 0, 0, - omega_oz, 0, np.sin(heading), -np.cos(heading)], # Phi_oy
                  #   [np.tan(phi)/(R+height), 0, 0, omega_oy, -omega_ox, 0, 0, 0], # Delta dot Phi_oz
                    [0, 0, 0, 0, 0, 0, 0], # Delta omega_x
                    [0, 0, 0, 0, 0, 0, 0]]) # Delta omega_y
                  #   [0, 0, 0, 0, 0, 0, 0, 0, -beta]]) # Delta omega_z
      # '''
      '''
      # модель из уравнений ошибок, упрощенная, с без курсового дрейфа, размерность 8
      A = np.array([[V_oy/(R+height)*np.tan(phi) - V_oz/(R+height), V_ox/(R+height)*np.tan(phi) + 2*U*np.sin(phi), -(V_ox/(R+height) + 2*U*np.cos(phi)), 0, -a_oz, a_oy, 0, 0], # Delta dot V_ox
                    [-2*(V_ox/(R+height)*np.tan(phi) + U*np.sin(phi)), -V_oz/(R+height), -V_oy/(R+height), a_oz, 0, -a_ox, 0, 0], # Delta dot V_oy
                    [2*(V_ox/(R+height) + U*np.cos(phi)), 2*V_oy/(R+height), 0, -a_oy, a_ox, 0, 0, 0], # Delta dot V_oz
                    [0, -1/(R+height), 0, 0, omega_oz, -omega_oy, 0, 0], # Phi_ox
                    [1/(R+height), 0, 0, - omega_oz, 0, omega_ox, 0, 0], # Phi_oy
                    [np.tan(phi)/(R+height), 0, 0, omega_oy, -omega_ox, 0, 0, 0], # Delta dot Phi_oz
                    [0, 0, 0, -beta, 0, 0, 0, 0], # Delta omega_x
                    [0, 0, 0, 0, -beta, 0, 0, 0]]) # Delta omega_y
                  #   [0, 0, 0, 0, 0, 0, 0, 0, -beta]]) # Delta omega_z
      # '''
      '''
      # модель из уравнений ошибок, упрощенная, с дрейфами, размерность 9
      A = np.array([[V_oy/(R+height)*np.tan(phi) - V_oz/(R+height), V_ox/(R+height)*np.tan(phi) + 2*U*np.sin(phi), -(V_ox/(R+height) + 2*U*np.cos(phi)), 0, -a_oz, a_oy, 0, 0, 0], # Delta dot V_ox
                    [-2*(V_ox/(R+height)*np.tan(phi) + U*np.sin(phi)), -V_oz/(R+height), -V_oy/(R+height), a_oz, 0, -a_ox, 0, 0, 0], # Delta dot V_oy
                    [2*(V_ox/(R+height) + U*np.cos(phi)), 2*V_oy/(R+height), 0, -a_oy, a_ox, 0, -0, 0, 0], # Delta dot V_oz
                    [0, -1/(R+height), 0, 0, omega_oz, -omega_oy, 1, 0, 0], # Phi_ox
                    [1/(R+height), 0, 0, - omega_oz, 0, omega_ox, 0, 1, 0], # Phi_oy
                    [np.tan(phi)/(R+height), 0, 0, omega_oy, -omega_ox, 0, 0, 0, 1], # Delta dot Phi_oz
                    [0, 0, 0, 0, 0, 0, -beta, 0, 0], # Delta omega_x
                    [0, 0, 0, 0, 0, 0, 0, -beta, 0], # Delta omega_y
                    [0, 0, 0, 0, 0, 0, 0, 0, -beta]]) # Delta omega_z
      # '''
      '''
      #модель из уравнений ошибок, расширенная дрейфами размерность 12
      A = np.array([[V_oy/(R+height)*np.tan(phi) - V_oz/(R+height), V_ox/(R+height)*np.tan(phi) + 2*U*np.sin(phi), -(V_ox/(R+height) + 2*U*np.cos(phi)), 0, -a_oz, a_oy, (V_ox/((R+height)*np.cos(phi)**2) + 2*U*np.cos(phi) )* V_oy +  2*U*np.sin(phi)*V_oz, 0, 0, 0, 0, 0], # Delta dot V_ox
                    [-2*(V_ox/(R+height)*np.tan(phi) + U*np.sin(phi)), -V_oz/(R+height), -V_oy/(R+height), a_oz, 0, -a_ox, -(V_ox/((R+height)*np.cos(phi)**2) + 2*U*np.cos(phi))*V_ox, 0, 0, 0, 0, 0], # Delta dot V_oy
                    [2*(V_ox/(R+height) + U*np.cos(phi)), 2*V_oy/(R+height), 0, -a_oy, a_ox, 0, -2*U*np.sin(phi)*V_ox, 0, -(k_2 - 2*nu**2), 0, 0, 0], # Delta dot V_oz
                    [0, -1/(R+height), 0, 0, omega_oz, - omega_oy, 0, 0, 0, 1, 0, 0], # Delta dot Phi_ox
                    [1/(R+height), 0, 0, -omega_oz, 0, omega_ox, -U*np.sin(phi), 0, 0, 0, 1, 0], # Delta dot Phi_oy
                    [np.tan(phi)/(R+height), 0, 0, omega_oy, -omega_ox, 0, U*np.cos(phi) + V_ox/((R+height)*np.cos(phi)**2), 0, 0, 0, 0, 1], # Delta dot Phi_oz
                    [0, 1/(R+height), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], # Delta dot varphi
                    [1/((R+height)*np.cos(phi)), 0, 0, 0, 0, 0, 0, 0, np.tan(phi)*V_ox/((R+height)*np.cos(phi)), 0, 0, 0], # Delta dot lambda
                    [0, 0, 1, 0, 0, 0, 0, 0, -k_1, 0, 0, 0], # \Delta dot h
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, -beta, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -beta, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -beta]])
      # '''
      filter.Phi = I + A*tau  # F - матрица процесса - размер dim_x на dim_x 
      filter.predict()
      filter.update((V_ins - V_o_gps).reshape(dim_z, 1))
      omega_est[:, i] = filter.x[-2:]
      x_est[:, i ] = filter.x
      str_for_write = str(data["Timestamp"][i]) + "," + str(float(filter.x[-2])) + "," + str(float(filter.x[-1])) + "\n"
      est_file.write(str_for_write)

'''
Сохранение оценок дрейфов (и нетолько) в файл
'''
vars = {'omega_est':omega_est}
with open("./Data_files/estimations.pkl", "wb") as file:
    pickle.dump(vars, file)

'''данные для оформления графиков'''
font_size = 20
line_width = 5

'''
строит оценки дрейфа
'''
# Скорости
# time = np.linspace(0, t_nav/60, N)
fig = go.Figure()
fig.add_trace(go.Scatter(x = time, y = x_est[0,:], name="$V_E$", line = dict(width = line_width)))
fig.add_trace(go.Scatter(x = time, y = x_est[1,:], name="$V_N$", line = dict(width = line_width)))
fig_name = "Оценки скоростей"
fig.update_layout(title= dict(text = fig_name, font=dict(size=font_size)) ,
                  xaxis_title="Время, мин",
                  yaxis_title="м/с",
                  legend = dict(orientation = "h", font = dict(size = font_size)))
# fig.write_image(f"./Графики для РПЗ и Листов/Моделирование_уравнений/Вариант_{num_variant}/Estimations_velocity_error_by_kalman_Variant_{num_variant}.png", width = 1598, height = 843)
fig.show()


# Углы
# time = np.linspace(0, t_nav/60, N)
fig = go.Figure()
fig.add_trace(go.Scatter(x = time, y = np.rad2deg(x_est[3,:])*3600, name="$\Phi_{x}$", line = dict(width = line_width)))
fig.add_trace(go.Scatter(x = time, y = np.rad2deg(x_est[4,:])*3600, name="$\Phi_{y}$", line = dict(width = line_width)))
# fig.add_trace(go.Scatter(x = time, y = x_est[-1,:], name="$\delta \omega_{oz}$"))
fig_name = "Оценки углов"
fig.update_layout(title= dict(text = fig_name, font=dict(size=font_size)) ,
                  xaxis_title="Время, мин",
                  yaxis_title="Угл. сек.",
                  legend = dict(orientation = "h", font = dict(size = font_size)))
# fig.write_image(f"./Графики для РПЗ и Листов/Моделирование_уравнений/Вариант_{num_variant}/Estimations_angle_palt_error_by_kalman_Variant_{num_variant}.png", width = 1598, height = 843)
fig.show()


# дрейф
# time = np.linspace(0, t_nav/60, N)
fig = go.Figure()
fig.add_trace(go.Scatter(x = time, y = np.rad2deg(x_est[-2,:])*3600, name="$\delta \omega_{bx}$", line = dict(width = line_width)))
fig.add_trace(go.Scatter(x = time, y = np.rad2deg(x_est[-1,:])*3600, name="$\delta \omega_{by}$", line = dict(width = line_width)))
# fig.add_trace(go.Scatter(x = time, y = x_est[-1,:], name="$\delta \omega_{oz}$"))
fig_name = "Оценки дрейфов"
fig.update_layout(title= dict(text = fig_name, font=dict(size=font_size)) ,
                  xaxis_title="Время, мин",
                  yaxis_title="Град/ч",
                  legend = dict(orientation = "h", font = dict(size = font_size)))
# fig.write_image(f"./Графики для РПЗ и Листов/Моделирование_уравнений/Вариант_{num_variant}/Estimations_bias_by_kalman_Variant_{num_variant}.png", width = 1598, height = 843)
fig.show()


print("Done!")
