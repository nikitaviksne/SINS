import numpy as np
import scipy
import matplotlib.pyplot as plt

'''
настройка отображения графиков
'''
line_width = 4
font_size = 30



g = 9.81
U = np.radians(15)/3600 #радианы в секунды
# U = 15
R_e = 6400*10**3 #м экваториальный радиус
R = R_e

nu = np.sqrt(g/R_e) #частота Шуллера
h = 1.0/100 #шаг соответствующий частоте 100 Гц

'''
Для скоростной коррекции
'''
print("Скоростная коррекция")
k1 = lambda w: 1.4*w
k2 = lambda w: w**2/nu**2 - 1

T_k_a = 1
beta_a = 1/T_k_a
T_k_g = 2
beta_w = 1/T_k_g
std_w = np.deg2rad(0.02)/3600#/np.sqrt(h) #радианы в секунду. 0.1 градус в час. Дисперсия гироскопов
std_a = 0.1*g*1e-3#/np.sqrt(h) #дисперсия акселерометров
std_V = 0.05/np.sqrt(h)

'''
# Мои значения в ДЗ
beta_w = 0.5
beta_a = 1
std_w = np.deg2rad(0.002)/3600 #радианы в секунду. 0.1 градус в час
std_a = 0.1*10**(-3)*g
std_V = 0.2
# '''

S_gyr = lambda w: 2 * beta_w * std_w**2/(np.pi*(w**2 + beta_w**2)) #спектральная плтность шума от дрейфа гироса
S_acc = lambda w: 2 * beta_a * std_a**2/(np.pi*(w**2 + beta_a**2)) #спектральная плтность шума от дрейфа акса
S_SNS_v = lambda w: std_V**2/(np.pi) #спектральная плтность шума от SNS

W_F_a_2 = lambda w, wx: (k2(wx)+1)**2/((R*nu**2 - R*w**2 + R*k2(wx)*nu**2)**2 + (R*w*k1(wx))**2) #квадрат передаточной функции от ошибки акса к ошибке вертикали
# W_F_om_2 = lambda w, wx: (w**2 + k1(wx)**2)/(R**2*(w**4+w**2*k1(w_x)**2 - 2*w**2*k2(w_x)*nu**2 - 2*w**2*nu**2 + k2(w_x)**2*nu**4 - 2*k2(w_x)*nu**4 + nu**4)) #квадрат передаточной функции от ошибки гиооса к ошибке вертикали
W_F_om_2 = lambda w, wx: (w**2 + k1(wx)**2)/((w**2 - (1 + k2(wx))*nu**2)**2 + w**2*k1(wx)**2) #квадрат передаточной функции от ошибки гиооса к ошибке вертикали
W_F_V_2 = lambda w, wx: (k1(wx)**2 + (w*k2(wx))**2)/((R*(k2(wx)+1)*nu**2 - R*w**2)**2 + (R*k1(wx)*w)**2) #квадрат передаточной функции от ошибки SNS к ошибке вертикали

S_w = lambda w, wx: W_F_om_2(w, wx)*S_gyr(w) #спектральная плотность ошибки вертикали от дрейфа гирсокопа
S_a = lambda w, wx: W_F_a_2(w, wx)*S_acc(w) #спектральная плотность ошибки вертикали от дрейфа акселерометра
S_V = lambda w, wx: W_F_V_2(w, wx)*S_SNS_v(w) #спектральная плотность ошибки вертикали от ошибки СНС

num_size = 300
sigma_sum = np.zeros(num_size)  #массив суммарной СКО
omega_st = 0 #начальная частота инегрирования
omega_fin =  100 #конечная частота инегрирования



D = np.zeros((4, num_size)) #массив для дисперсий по дрейфу гироскопа (w), акселерометра (a), координате (phi) и сумма дисперсий

'''
Эту оптимальную частоту надо поставлять в k1,k2,k3, получать числа, эти численные значения
подставлять в спектральные плотности, спектральная плотность, тогда получится функцией только от
частоты, и чтобы найти дисперсию, надо интегрировать эти спектральные плотности по единственной
переменной -- частоте
'''

"""Когда спектральная плотность - функция только одной переменной, омега, она же подставляется в
коэффициенты коррекции. А интегрирование не производится. Тут вопрос, почему по теории надо интегрировать,
а на практике не надо.

И еще, совершенно точно надо переводить из радиан в градусы под квадратом, т.е *(rad2deg*3600)**2
"""


# Первый вариант (интегрирование другим методом)
num_samples = int((omega_fin - omega_st)/h)
x_for_integr = np.linspace(omega_st, omega_fin, num_samples)


x = np.linspace(1*nu, 4*nu, num_size)
for num, w_x in enumerate(x):
    D[0, num] = scipy.integrate.trapezoid(S_w(x_for_integr, w_x), dx=0.01)*(np.degrees(1)*3600)**(2)
    D[1, num] = scipy.integrate.trapezoid(S_a(x_for_integr, w_x), dx=0.01)*(np.degrees(1)*3600)**(2)
    D[2, num] = scipy.integrate.trapezoid(S_V(x_for_integr, w_x), dx=0.01)*(np.degrees(1)*3600)**(2)
    D[3, num] = (D[0][num] + D[1][num] + D[2][num])
    sigma_sum[num] = np.sqrt(D[3][num])

min_std_Phi = min(sigma_sum)
print('min_std_phi = ', min_std_Phi) #Сам минимальный элемент массива суммарного СКО
index_min = np.argmin(sigma_sum)
freq_min_std_Phi = x[index_min]#оптимальная собственная частота, соответсвтующая минимальному суммарному СКО
print("Индекс минимального элемента = ", index_min)
print("Оптимальная собственная частота = ", freq_min_std_Phi)
print("Оптимальная относительная частота = ", freq_min_std_Phi/nu)
#рассчитываем коэффициенты k1, k2, (k3 если есть)
k1_optima = k1(freq_min_std_Phi)
k2_optima = k2(freq_min_std_Phi)
print('k1_optima = ', k1_optima)
print('k2_optima = ', k2_optima)


x_plot = x/nu

fig1, (fig1_ax1) = plt.subplots(1,1)
fig1.suptitle("Дисперсия")
fig1_ax1.plot(x_plot, D[0], label = "Дисперсия от дрефа гироскопа", linewidth = line_width);
fig1_ax1.plot(x_plot, D[1], label = "Дисперсия от дрефа акселерометра", linewidth = line_width);
fig1_ax1.plot(x_plot, D[2], label = "Дисперсия от ошибки СНС по скорости", linewidth = line_width);
fig1_ax1.legend(loc = 'best');
fig1_ax1.grid(True);



fig2, (fig2_ax1) = plt.subplots(1,1)
fig2.suptitle("Суммарное СКО")
fig2_ax1.plot(x_plot, sigma_sum, label = "Суммарное СКО", linewidth = line_width);
fig2_ax1.legend(loc = 'best');
fig2_ax1.grid(True);

plt.show()
exit();
'''
Для позиционной коррекции
'''
print("Позиционная коррекция")

k1 = lambda w: 1.75*w
k2 = lambda w: 2.15*w**2/nu**2 - 1
k3 = lambda w: w**3/nu**2 - 1.75*w

# пример

# beta_w = 0.1
# beta_a = 0.2
# std_w = np.deg2rad(0.1)/3600#/np.sqrt(h) #радианы в секунду. 0.1 градус в час
# std_a = 5e-5*g#/np.sqrt(h)
std_phi = 3/R

# # Мои значения в ДЗ
# beta_w = 1
# beta_a = 0.5
# std_w = np.radians(0.001)/3600#/np.sqrt(h) #радианы в секунду. 0.01 градус в час
# std_a = 0.1*10**(-5)*g#/np.sqrt(h)
# std_phi = 3/R#/np.sqrt(h) #0.1 Метра, СКО ошибки по координате от СНС

S_gyr = lambda w: 2*beta_w*std_w**2/(np.pi*(w**2+beta_w**2)) #спектральная плтность шума от дрейфа гироса
S_acc = lambda w: 2*beta_a*std_a**2/(np.pi*(w**2+beta_a**2)) #спектральная плтность шума от дрейфа акса
S_SNS_phi = lambda w: std_phi**2/(np.pi) #спектральная плтность шума от SNS

W_F_a_2 = lambda w, wx: (w**2 + (k1(wx) + k3(wx))**2 )/(( ( (1 + k2(wx) )*nu**2*w - w**3 )**2 + ( (k1(wx) + k3(wx) )*nu**2 - k1(wx)*w**2)**2)*R**2) #квадрат передаточной функции от ошибки акса к ошибке вертикали
W_F_om_2 = lambda w, wx: ((k2(wx)*nu**2 - w**2)**2 + k1(wx)**2*w**2)/( ( (1 + k2(wx) )*nu**2*w - w**3 )**2 + ( (k1(wx) + k3(wx) )*nu**2 - k1(wx)*w**2)**2) #квадрат передаточной функции от ошибки гиооса к ошибке вертикали
W_F_phi_2 = lambda w, wx: (k3(wx)**2 * w**4 + k2(wx)**2*nu**4*w**2)/( ( (1 + k2(wx) )*nu**2*w - w**3 )**2 + ( (k1(wx) + k3(wx) )*nu**2 - k1(wx)*w**2)**2) #квадрат передаточной функции от ошибки SNS к ошибке вертикали

S_w = lambda w, w_x: W_F_om_2(w, w_x)*S_gyr(w) #спектральная плотность ошибки вертикали от дрейфа гирсокопа
S_a = lambda w, w_x: W_F_a_2(w, w_x)*S_acc(w) #спектральная плотность ошибки вертикали от дрейфа акселерометра
S_phi = lambda w, w_x: W_F_phi_2(w, w_x)*S_SNS_phi(w) #спектральная плотность ошибки вертикали от ошибки СНС

num_size = 500 #количество точек при выводе графика

D = np.zeros((4, num_size))#массив для дисперсий по дрейфу гироскопа (w), акселерометра (a), координате (phi) и сумма дисперсий
sigma_sum = np.zeros(num_size)  #массив суммарной СКО

omega_st = 0 #начальная частота инегрирования
omega_fin =  200 #конечная частота инегрирования
# omega_fin =  np.inf #конечная частота инегрирования

# интегрирование трапециями
num_samples = int((omega_fin - omega_st)/h) #количество точек для интегрирования
x_for_integr = np.linspace(omega_st, omega_fin, num_samples)

x = np.linspace(3*nu, 20*nu,  num_size)

for num, w_x in enumerate(x):
    D[0, num] = scipy.integrate.trapezoid(S_w(x_for_integr, w_x), x_for_integr)*(np.degrees(1)*3600)**(2)
    D[1, num] = scipy.integrate.trapezoid(S_a(x_for_integr, w_x), x_for_integr)*(np.degrees(1)*3600)**(2)
    D[2, num] = scipy.integrate.trapezoid(S_phi(x_for_integr, w_x), x_for_integr)*(np.degrees(1)*3600)**(2)
    D[3, num] = (D[0][num] + D[1][num] + D[2][num])
    sigma_sum[num] = np.sqrt(D[3][num])

min_std_Phi = min(sigma_sum)
print('min_std_phi = ', min_std_Phi) #Сам минимальный элемент массива суммарного СКО
index_min = np.argmin(sigma_sum)
freq_min_std_Phi = x[index_min]#собственная частота, соответсвтующая минимальному суммарному СКО
print("index_min = ", index_min)
print("min_freq = ",freq_min_std_Phi)
print("Оптимальная относительная частота = ", freq_min_std_Phi/nu)
#рассчитываем коэффициенты k1, k2, (k3 если есть)
k1_optima = k1(freq_min_std_Phi)
k2_optima = k2(freq_min_std_Phi)
k3_optima = k3(freq_min_std_Phi)
print('k1_optima = ', k1_optima)
print('k2_optima = ', k2_optima)
print('k3_optima = ', k3_optima)

x_plot = x/nu

fig = go.Figure()
fig.add_trace(go.Scatter(x = x_plot, y = D[0], name="Дисперсия от дрефа гироскопа", line = dict(width = line_width)))
fig.add_trace(go.Scatter(x = x_plot, y = D[1], name="Дисперсия от дрефа акселерометра", line = dict(width = line_width)))
fig.add_trace(go.Scatter(x = x_plot, y = D[2], name="Дисперсия от ошибки СНС по координате", line = dict(width = line_width)))
fig.add_trace(go.Scatter(x = x_plot, y = D[3], name="Суммарная дисперсия", line = dict(width = line_width*3)))
fig.update_layout(title= dict(text = "Дисперсия ошибки ориентации", font = dict(size=font_size)),
                  xaxis=dict(title="Относительная частота", titlefont=dict(size=font_size), tickfont = dict(size = font_size)),
                  yaxis=dict(title="угл. с. ^2", titlefont=dict(size=font_size), tickfont = dict(size = font_size)),
                  legend=dict(
    orientation="h",
    font=dict(size=font_size-2)
    # yanchor="bottom",
    # entrywidth=300,
    # y=10,
    # xanchor="right",
    # x=-0 # Смещение начальной точки в направлении отсчета xanchor
))
fig.show()


fig = go.Figure()
fig.add_trace(go.Scatter(x = x_plot, y = sigma_sum, name = "Суммарное СКО", line = dict(width = line_width)))
fig.update_layout(title = dict(text = "Суммарное СКО ошибки ориентации", font = dict(size=font_size)),
                  xaxis=dict(title="Относительная частота", titlefont=dict(size=font_size), tickfont = dict(size = font_size)),
                  yaxis=dict(title="угл. с.", titlefont=dict(size=font_size), tickfont = dict(size = font_size)))
fig.show()

exit(0)

'''Вертикальный канал'''
print("Вертикальный канал")
k1 = lambda w: 1.4*w
k2 = lambda w: w**2 + 2*nu**2

# beta_w = 0.1
# beta_a = 5
# std_w = np.radians(0.1)/3600#/np.sqrt(h) #радианы в секунду. 0.1 градус в час
# std_a = 1e-3*g#/np.sqrt(h)
# mu_a = 1e-3*g # постоянная сосотавляющая дрейфа акселерометра
# std_V = 0.05#/np.sqrt(h)

# Мои значения в ДЗ
beta_w = 0.5
beta_a = 1
std_a = 1*g*10**(-3)# 1e-3*g
std_h = 3
mu_a = 1e-3*g

S_acc = lambda w: 2 * beta_a * std_a**2/(np.pi*(w**2 + beta_a**2)) #спектральная плтность шума от дрейфа акса
S_dh = lambda w: std_h**2/(np.pi) #спектральная плтность шума от корректора по высоте

W_H_a_2 = lambda w, wx: (1)/((-w**2 +k2(wx) - 2*nu**2)**2 + (w*k1(wx)**2)) #квадрат передаточной функции от ошибки акса к ошибке вертикали
# W_F_om_2 = lambda w, wx: (w**2 + k1(wx)**2)/(R**2*(w**4+w**2*k1(w_x)**2 - 2*w**2*k2(w_x)*nu**2 - 2*w**2*nu**2 + k2(w_x)**2*nu**4 - 2*k2(w_x)*nu**4 + nu**4)) #квадрат передаточной функции от ошибки гиооса к ошибке вертикали
W_H_h_2 = lambda w, wx: (k2(wx)**2 + w**2*k1(wx)**2)/(w**4 + (-2*k2(wx) + 4*nu**2 + k1(wx)**2)*w**2 + (k2(wx) - 2*nu**2)**2) #квадрат передаточной функции от ошибки SNS к ошибке вертикали

S_mu_a = lambda w: mu_a/w**2
S_a = lambda w, wx: W_F_a_2(w, wx)*S_acc(w) #спектральная плотность ошибки вертикали от дрейфа акселерометра
S_h = lambda w, wx: W_H_h_2(w, wx)*S_dh(w) #спектральная плотность ошибки вертикали от ошибки СНС

num_size = 300
sigma_sum = np.zeros(num_size)  #массив суммарной СКО
omega_st = 0 #начальная частота инегрирования
omega_fin =  100 #конечная частота инегрирования



D = np.zeros((4, num_size)) #массив для дисперсий по дрейфу гироскопа (w), акселерометра (a), координате (phi) и сумма дисперсий

num_samples = int((omega_fin - omega_st)/h)
x_for_integr = np.linspace(omega_st, omega_fin, num_samples)


x = np.linspace(80*nu, 150*nu, num_size)
for num, w_x in enumerate(x):
    D[0, num] = S_mu_a(w_x)**2
    D[1, num] = scipy.integrate.trapezoid(S_a(x_for_integr, w_x), x_for_integr)#*(np.degrees(1)*3600)**(2)
    D[2, num] = scipy.integrate.trapezoid(S_h(x_for_integr, w_x), x_for_integr)#*(np.degrees(1)*3600)**(2)
    D[3, num] = (D[0][num] + D[1][num] + D[2][num])
    sigma_sum[num] = np.sqrt(D[3][num])

min_std_h = min(sigma_sum)
print('min_std_h = ', min_std_h) #Сам минимальный элемент массива суммарного СКО
index_min = np.argmin(sigma_sum)
freq_min_std_Phi = x[index_min]#оптимальная собственная частота, соответсвтующая минимальному суммарному СКО
print("Индекс минимального элемента = ", index_min)
print("Оптимальная собственная частота = ", freq_min_std_Phi)
print("Оптимальная относительная частота = ", freq_min_std_Phi/nu)
#рассчитываем коэффициенты k1, k2, (k3 если есть)
k1_optima = k1(freq_min_std_Phi)
k2_optima = k2(freq_min_std_Phi)
print('k1_optima = ', k1_optima)
print('k2_optima = ', k2_optima)


x_plot = x/nu

fig = go.Figure()
fig.add_trace(go.Scatter(x = x_plot, y = D[0]-0*6.28e-6, name="Дисперсия постоянной составляющей акселерометра", line = dict(width = line_width)))
fig.add_trace(go.Scatter(x = x_plot, y = D[1], name="Дисперсия от дрефа акселерометра", line = dict(width = line_width)))
fig.add_trace(go.Scatter(x = x_plot, y = D[2]-0*6.28e-6, name="Дисперсия от ошибки корректора по высоте", line = dict(width = line_width)))
fig.update_layout(title= dict(text = "Дисперсия ошибки высоты", font = dict(size=font_size)),
                  xaxis=dict(title="Относительная частота", titlefont=dict(size=font_size), tickfont = dict(size = font_size)),
                  yaxis=dict(title="м^2", titlefont=dict(size=font_size), tickfont = dict(size = font_size)),
                  legend=dict(
                orientation="h",
                font=dict(size=font_size-2)
                # yanchor="bottom",
                # entrywidth = 300,
                # y=10,
                # xanchor="right",
                # x=-0 # Смещение начальной точки в направлении отсчета xanchor
                ))
fig.show()


fig = go.Figure()
fig.add_trace(go.Scatter(x = x_plot, y = sigma_sum-0*1.25, name = "Суммарное СКО", line = dict(width = line_width)))
fig.update_layout(title = dict(text ="Суммарное СКО ошибки высоты", font = dict(size=font_size)),
                  xaxis=dict(title="Относительная частота", titlefont=dict(size=font_size), tickfont = dict(size = font_size)),
                  yaxis=dict(title="м", titlefont=dict(size=font_size), tickfont = dict(size = font_size)))
fig.show()
