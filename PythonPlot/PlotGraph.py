import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import LineCollection
import sys #для приема имени файла для построения через аргумент команды
import os
import logging
import argparse

def set_size(width_pt, fraction=1):
    """
    Рассчитывает figsize в дюймах.
    width_pt: значение \textwidth из LaTeX
    fraction: доля ширины (например, 0.9)
    """
    fig_width_pt = width_pt * fraction
    inches_per_pt = 1 / 72.27
    
    # Золотое сечение для высоты
    golden_ratio = (5**.5 - 1) / 2
    
    fig_width_in = fig_width_pt * inches_per_pt
    fig_height_in = fig_width_in * golden_ratio
    
    return (fig_width_in, fig_height_in)

def saveFigure(fig, full_path, format_type):
    """Сохраняет фигуру. Директория должна существовать."""
    try:
        os.makedirs(os.path.dirname(full_path), exist_ok=True)
        fig.savefig(full_path+f".{format_type}", bbox_inches='tight', pad_inches=0.01, format=format_type)
        logging.info(f"График сохранен: {full_path}")
    except Exception as e:
        logging.error(f"Не удалось сохранить график {full_path}: {e}")
        
def apply_lowpass_filter(data_series, time_series, tau=1e6, alpha=1e-5, dt=None):
    """
    Применяет фильтр нижних частот к данным.
    Возвращает отфильтрованный сигнал и его скользящее среднее.
    """
    if dt is None:
        # Вычисляем средний шаг по времени, если не передан
        dt = np.mean(np.diff(time_series))
    
    flf = [data_series.iloc[0]]  # Начальное значение фильтра
    flf_alpha = [data_series.iloc[0]]
    mean = [data_series.iloc[0]]

    for i in range(1, len(data_series)):
        x = data_series.iloc[i]
        # Рекурсивный фильтр с постоянной времени tau
        flf.append((x + tau/dt * flf[-1]) / (1 + tau/dt))
        # Экспоненциальное сглаживание с коэффициентом alpha
        flf_alpha.append(flf_alpha[-1] * (1 - alpha) + x * alpha)
        # Скользящее среднее
        mean.append(mean[-1] + (flf[-1] - mean[-1]) / i)
        
    return flf, flf_alpha, mean

plt.rcParams.update({
    'figure.dpi': 150,       # Экранное разрешение
    'savefig.dpi': 600,      # Разрешение при сохранении выбор значения DPI зависит от назначения графика: 
# 72–96 DPI — для веб-страниц и экранного отображения; 
# 150 DPI — для базовых презентаций и внутренних документов; 
# 300 DPI — для публикаций, профессиональных отчётов и печати; 
# 600+ DPI — для научных публикаций и высококачественной печати. 
    'font.size': 9,
    'font.family': 'serif',
    'axes.titlesize': 8,
    'axes.labelsize': 7,
    'legend.fontsize': 7,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'lines.linewidth': 1.0,
    'lines.markersize': 4,
    'figure.titlesize': 11,   # suptitle крупнее
    'axes.grid': True,       # Включаем сетку по умолчанию
    'grid.alpha': 0.6,
    'grid.linestyle': '--',
    'legend.loc': 'best',
    'lines.linewidth': 1.5,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02
})

logging.basicConfig(level = logging.CRITICAL)
# logging.basicConfig(level = logging.DEBUG)

parser = argparse.ArgumentParser(description='Plotting Graphs')
parser.add_argument('NavRes', type=str, help='File of result INS modelling')
parser.add_argument('KalmanEst', type=str, help='File of results Kalman estimation')
parser.add_argument('--SaveDir', type=str, help='Directory to save the figure (optional). If not will not saving')
parser.add_argument('--PGF', action='store_true', help='Enable PGF backend for LaTeX')
# Теперь использование: script.py data.csv est.csv --PGF
# args.PGF будет True, если флаг передан, и False в противном случае
args = parser.parse_args()

logging.log(logging.DEBUG, f"Agument SaveDir is {args.SaveDir}")
logging.log(logging.DEBUG, f"Agument PGF is {args.PGF}")


try:
    if args.PGF:# 1. Настройка бэкенда для PGF
        mpl.use("pgf")
except RuntimeError:
    logging.log(logging.DEBUG, "Вероятно, не установлен LuaLaTeX в системе")
try:

    if args.PGF:# 2. Установка параметров текста (чтобы шрифты совпадали с LuaLaTeX)
        mpl.rcParams.update({
            "pgf.texsystem": "lualatex",
            "font.family": "serif",  # или 'sans-serif'
            "text.usetex": True,
            "pgf.rcfonts": False,    # не использовать шрифты matplotlib
        })
except RuntimeError:
    logging.log(logging.DEBUG, "Вероятно, не установлен LuaLaTeX в системе")

width = 500.484 # размер \showthe\textwidth из Latex 500.484pt.

logging.log(logging.DEBUG, f"Size in python = {set_size(width)}")

path_nav_sol = sys.argv[1]
path_est = sys.argv[2]

savefig = False # не сохраняем графики
if (len(sys.argv) > 3): #если не указал третим параметром путь сохранения графиков, то оставляем переменную savefig в значении false
    savefig = True
    savefig_path_dir = args.SaveDir # путь к дериктории
# data =pd.read_csv("~/InertialNavigation/data/"+relative_path, delimiter=";");
data = pd.read_csv(path_nav_sol, delimiter=";");
estimations = pd.read_csv(path_est, delimiter=",");
# После загрузки data и estimations
required_data_cols = ["d_VE", 
                      "d_VN", 
                      "d_VUp", 
                      "Heading", 
                      "Roll", 
                      "Pitch",
                      ]
for col in required_data_cols:
    if col not in data.columns:
        logging.error(f"В файле навигационного решения отсутствует столбец: {col}")
        # Можно либо выйти, либо создать столбец с NaN, чтобы график не падал
        data[col] = np.nan
start = 0;
stop = -1; #1*60*100

time = np.linspace(0, (data.iloc[:, 0].size - 1)/100/60, data.iloc[:, 0].size) #в минутах /100/60
# time = np.linspace(0, (data.iloc[:, 0].size - 1), data.iloc[:, 0].size) # в тактах

round = 100

'''Ошибки горизонтальных линейный скоростей'''
# по совету DeepSeek
# Создаем фигуру 2x3
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=set_size(width, 1.2), sharex=True) 
# sharex=True очень полезно для временных рядов

# --- Ряд 1: Ошибки скорости ---
# Восточная скорость
ax = axes[0, 0]
ax.plot(time, data["d_VE"], label="Ошибка БИНС", linewidth=1.0)
ax.plot(time, estimations["Ve"], '--', label="Оценка ФК", linewidth=1.5)
ax.set_ylabel("Ve, м/с")
ax.set_title("Восточная скорость")
ax.legend(fontsize=9)

# Северная скорость
ax = axes[0, 1]
ax.plot(time, data["d_VN"], label="Ошибка БИНС")
ax.plot(time, estimations["Vn"], '--', label="Оценка ФК")
ax.set_ylabel("Vn, м/с")
ax.set_title("Северная скорость")
ax.legend(fontsize=9)

# --- Ряд 2: Скорости дрейфов гироскопов ---
# Перед этим блоком вычислите flf_alpha_x, flf_alpha_y, flf_alpha_z
flf_alpha_x, _, _ = apply_lowpass_filter(estimations["d_omega_x"], time, dt=100)
flf_alpha_y, _, _ = apply_lowpass_filter(estimations["d_omega_y"], time, dt=100)

# Канал X
ax = axes[1, 0]
ax.plot(time, np.rad2deg(estimations["d_omega_x"]) * 3600, alpha=0.7, label="Оценка")
ax.plot(time, np.rad2deg(flf_alpha_x) * 3600, linewidth=2, label="Сглаж. оценка")
ax.set_ylabel("град/час")
ax.set_xlabel("Время, мин")
ax.set_title("Дрейф гироскопа X")
ax.axhline(y=0.1, color='r', linestyle=':', alpha=0.5) # Значение модели для сравнения
ax.legend(fontsize=9)

# Канал Y
ax = axes[1, 1]
ax.plot(time, np.rad2deg(estimations["d_omega_y"]) * 3600, alpha=0.7, label = "Оценка")
ax.plot(time, np.rad2deg(flf_alpha_y) * 3600, linewidth=2, label="Сглаж. оценка")
ax.set_xlabel("Время, мин")
ax.set_title("Дрейф гироскопа Y")
ax.axhline(y=0.1, color='r', linestyle=':', alpha=0.5)
ax.legend(fontsize=9)


# Общий заголовок
fig.suptitle("Сравнение ошибок БИНС и оценок фильтра Калмана", fontsize=16, y=1.02)
# Настройка расстояний между subplots
fig.tight_layout()
if args.SaveDir:
    fig1_path = os.path.join(args.SaveDir, "Ошибки_горизонтальных_скоростей")
    saveFigure(fig, fig1_path, 'pdf' if not args.PGF else 'pgf')
# разворот на весь экран
# figManager = plt.get_current_fig_manager()
# figManager.window.showMaximized()
# fig1.canvas.manager.full_screen_toggle() # делаем полноэкранный режим, чтобы сохранялись кортинки в нормальном размере
# manager = plt.get_current_fig_manager()
# manager.window.state('zoomed')

'''Ошибки вертикальной линейной скорости'''
# по совету DeepSeek
# Создаем фигуру 2x3
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=set_size(width, 1.2), sharex=True) 
# sharex=True очень полезно для временных рядов

# --- Ряд 1: Ошибки скорости ---
# Вертикальная скорость
ax = axes[0]
ax.plot(time, data["d_VUp"], label="Ошибка БИНС")
# ax.plot(time, estimations["Vup"], '--', label="Оценка ФК") # Раскомментируйте, если есть
ax.set_ylabel("Vup, м/с")
ax.set_title("Вертикальная скорость")
ax.legend(fontsize=9)

# --- Ряд 2: Скорости дрейфов гироскопов ---
# Перед этим блоком вычислите flf_alpha_x, flf_alpha_y, flf_alpha_z
flf_alpha_z, _, _ = apply_lowpass_filter(estimations["d_omega_z"], time, dt=100)

# Канал Z
ax = axes[1]
ax.plot(time, np.rad2deg(estimations["d_omega_z"]) * 3600, alpha=0.7, label="Оценка")
ax.plot(time, np.rad2deg(flf_alpha_z) * 3600, linewidth=2, label="Сглаж. оценка")
ax.set_xlabel("Время, мин")
ax.set_title("Дрейф гироскопа Z")
# ax.axhline(y=..., color='r', linestyle=':')
ax.legend(fontsize=9)

# Общий заголовок
fig.suptitle("Сравнение ошибок БИНС и оценок фильтра Калмана", fontsize=16, y=1.02)
# Настройка расстояний между subplots
fig.tight_layout()
if args.SaveDir:
    fig1_path = os.path.join(args.SaveDir, "Ошибка_вертикальной_скорости")
    saveFigure(fig, fig1_path, 'pdf' if not args.PGF else 'pgf')

'''Углы ориентации'''
fig2, (fig2_ax1, fig2_ax2, fig2_ax3) = plt.subplots(nrows = 1, ncols = 3, 
                                                    figsize=set_size(width, fraction=0.9),
                                                    sharey='row',
                                                    )
fig2.suptitle("Углы")
# Угол курса
fig2_ax1.set_title("Угол курса");
fig2_ax1.set_xlabel("мин")
fig2_ax1.set_ylabel("Угл. мин.")
fig2_ax1.plot(time, np.round(np.rad2deg(data["Heading"])*60,round), label="Курс");
fig2_ax1.plot(time, np.round(np.rad2deg(estimations["d_Psi"]),round), label="$\hat{\Delta\Psi}$");
fig2_ax1.legend(loc="lower center")
fig2_ax1.grid(True)
# Угол крена
fig2_ax2.set_title("Угол крена");
fig2_ax2.set_xlabel("Мин")
fig2_ax2.set_ylabel("Угл. мин")
fig2_ax2.plot(time, np.round(np.rad2deg(data["Roll"])*60,round), label="$\gamma$");
fig2_ax2.plot(time, np.round(np.rad2deg(estimations["d_Roll"])*60,round), label="$\hat{\Delta\gamma}$");
fig2_ax2.set_ylim(-1.5, 2)
fig2_ax2.legend(loc="lower center")
fig2_ax2.grid(True)
# Угол тангажа
fig2_ax3.set_title("Угол тангажа");
fig2_ax3.set_xlabel("Мин")
fig2_ax3.set_ylabel("Угл. мин")
fig2_ax3.plot(time, np.round(np.rad2deg(data["Pitch"])*60,round), label="$\\theta$");
fig2_ax3.plot(time, np.round(np.rad2deg(estimations["d_Pitch"])*60,round), label="$\hat{\Delta\\theta}$");
fig2_ax3.set_ylim(-2, 1.5)
fig2_ax3.legend(loc="lower center")
fig2_ax3.grid(True)
if args.SaveDir:
    fig2_path = os.path.join(args.SaveDir, "Углы_ориентации")
    saveFigure(fig2, fig2_path, format_type='pdf' if not args.PGF else "pgf")

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
fig3, (fig3_ax1, fig3_ax2, fig3_ax3) = plt.subplots(nrows = 1, ncols = 3, 
                                                    figsize=set_size(width, fraction=0.9),
                                                    constrained_layout=True, # Более продвинутый алгоритм распределения места
                                                    sharey='row',
                                                    )
# fig3, (fig3_ax1, fig3_ax2, fig3_ax3) = plt.subplots(nrows = 1, ncols = 3, 
#                                                     figsize=set_size(1080, 1920), constrained_layout=True # Более продвинутый алгоритм распределения места
#                                                     )
# Восточное направление
fig3.suptitle("Оценка скоростей дрейфов гироскопов")
fig3_ax1.set_title("Оценка скоростей\nдрейфа гироскопа\nканала X")
fig3_ax1.set_xlabel("Мин")
fig3_ax1.set_ylabel("Град/час")
fig3_ax1.set_ylim(-0.25, 0.25)
'''
fig3_ax1.axhline(y = 0.05, xmin=0, color = "#FF0000", label = "Модель, body")
fig3_ax1.axhline(y = 0.05*np.cos(np.deg2rad(50)) + 0.05*np.sin(np.deg2rad(50)), xmin=0, color = "#db0db9", linestyle = "-.", label = "Модель, psi = 50$\circ$")
fig3_ax1.axhline(y = -0.05*np.sin(np.deg2rad(50)) + 0.05*np.cos(np.deg2rad(50)), xmin=0, color = "#db0db9",  linestyle = "-.", label = "Модель, psi = 50$\circ$")
fig3_ax1.axhline(y = 0.05*np.cos(np.deg2rad(50+90)) + 0.05*np.sin(np.deg2rad(50+90)), xmin=0, color = "#0ddb60", linestyle = "--", label = "Модель, psi = 140$\circ$")
fig3_ax1.axhline(y = -0.05*np.sin(np.deg2rad(50+90)) + 0.05*np.cos(np.deg2rad(50+90)), xmin=0, color = "#0ddb60",  linestyle = "--", label = "Модель, psi = 140$\circ$")
'''
fig3_ax1.plot(time, np.round(np.rad2deg(estimations["d_omega_x"]), round)*3600, label="$\delta\omega_x$");
fig3.set_size_inches(set_size(width))
# fig3_ax1.plot(time, np.round(np.rad2deg(mean[0]), round)*3600, label="$E(\omega_x^F)$");
# fig3_ax1.plot(time, np.round(np.rad2deg(flf[0]), round)*3600, label="$\omega_x^F$", alpha = 0.5, linestyle = "--", marker = "*");
fig3_ax1.plot(time, np.round(np.rad2deg(flf_alpha[0]), round)*3600, label="$\omega_{x}^{F\\alpha}$");
fig3_ax1.grid(True)
fig3_ax1.legend(loc="lower center")
#дрейф северного гироскопа
fig3_ax2.set_title("Оценка скоростей\nдрейфа гироскопа\nканала Y")
fig3_ax2.set_xlabel("Мин")
# fig3_ax2.set_ylabel("Град/час")
fig3_ax2.set_ylim(-0.25, 0.25)
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
fig3_ax2.legend(loc="lower center",)
#дрейф курсового гироскопа
# '''
fig3_ax3.set_title("Оценка скоростей\nдрейфа гироскопа\nканала Z")
fig3_ax3.set_xlabel("Мин")
# fig3_ax3.set_ylabel("Град/час")
fig3_ax3.set_ylim(-0.25, 0.25)
fig3_ax3.plot(time, np.round(np.rad2deg(estimations["d_omega_z"]), round)*3600, label="$\delta\omega_z$");
fig3_ax3.grid(True)
fig3_ax3.legend(loc="lower center")
# '''
if args.SaveDir:
    fig3_path = os.path.join(args.SaveDir, "Оценки_дрейфов_гироскопов")
    saveFigure(fig3, fig3_path, format_type='pdf' if not args.PGF else "pgf")

fig20, (fig20_ax1) = plt.subplots(ncols=1, nrows=1)
fig20.suptitle("Скорости ИНС")
fig20_ax1.plot(time, np.round(data["Ve"], round), label = "Ve")
fig20_ax1.plot(time, np.round(data["Vn"], round), label = "Vn")
fig20_ax1.set_xlabel("Мин")
fig20_ax1.set_ylabel("м/с")
fig20_ax1.grid(True)
fig20_ax1.legend(loc="lower center")

'''Малые приращения линейных скоростей'''
fig25, (fig25_ax1, fig25_ax2) = plt.subplots(1, 2)
# Восточная составляющая
fig25.suptitle("Малые приращения скоростей")
fig25_ax1.set_title("Малые приращения восточной скорости");
fig25_ax1.set_xlabel("мин")
fig25_ax1.set_ylabel("м/с")
fig25_ax1.plot(time, np.round(data["Woe"], round), linestyle="None", marker="+", label="Woe");
fig25_ax1.legend(loc="best")
fig25_ax1.grid(True)
# Северная составляющая
fig25_ax2.set_title("Малые приращения северной скорости");
fig25_ax2.set_xlabel("мин")
fig25_ax2.set_ylabel("м/с")
fig25_ax2.plot(time, np.round(data["Won"],round), linestyle="None", marker="+", label="Won");
fig25_ax2.legend(loc="best")
fig25_ax2.grid(True)


'''Траектория'''
fig35, ax = plt.subplots(figsize=set_size(width, 0.7))
fig35.suptitle("Траектория движения, окрашенная по времени")

lon = np.rad2deg(data["Lambda"])
lat = np.rad2deg(data["Phi"])

points = np.array([lon, lat]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# Нормализуем время для цвета
norm = plt.Normalize(time.min(), time.max())
lc = LineCollection(segments, cmap='viridis', norm=norm, linewidth=2, alpha=0.8)
lc.set_array(time)
line = ax.add_collection(lc)
ax.set_xlim(lon.min(), lon.max())
ax.set_ylim(lat.min(), lat.max())

# Добавляем colorbar
cbar = fig.colorbar(line, ax=ax, label='Время, мин', shrink=0.8)
# Подписываем начальную и конечную точки
ax.annotate('Старт', xy=(lon.iloc[0], lat.iloc[0]), xytext=(5, 5), textcoords='offset points', color='blue')
ax.annotate('Финиш', xy=(lon.iloc[-1], lat.iloc[-1]), xytext=(5, -15), textcoords='offset points', color='red')

ax.set_xlabel("Долгота, град")
ax.set_ylabel("Широта, град")
ax.grid(True, linestyle='--', alpha=0.7)
ax.set_aspect('equal', adjustable='box') # Очень важно для правильного отображения геоданных!
if args.SaveDir:
    fig35_path = os.path.join(args.SaveDir, "Траектория_с_цветом")
    saveFigure(fig35, fig35_path, format_type='pdf' if not args.PGF else "pgf")

'''Дрейфы акселерометров'''
fig40, (fig40_ax1, fig40_ax2, fig40_ax3) = plt.subplots(nrows = 1, ncols = 3, 
                                                    figsize=set_size(width, fraction=0.9),
                                                    constrained_layout=True, # Более продвинутый алгоритм распределения места
                                                    sharey='row',
                                                    )
# fig40, (fig40_ax1, fig40_ax2, fig40_ax3) = plt.subplots(nrows = 1, ncols = 3, 
#                                                     figsize=set_size(1080, 1920), constrained_layout=True # Более продвинутый алгоритм распределения места
#                                                     )
# Восточное направление
fig40.suptitle("Оценка скоростей дрейфов акселерометров")
fig40_ax1.set_title("Оценка скоростей\nдрейфа акселерометра \nканала X")
fig40_ax1.set_xlabel("Мин")
fig40_ax1.set_ylabel("м/(с^2)")
fig40_ax1.set_ylim(-1e-4, 1e-2)
fig40_ax1.plot(time, np.round(estimations["d_a_x"], round), label="$\delta a_x$");
fig40.set_size_inches(set_size(width))
# fig40_ax1.plot(time, np.round(np.rad2deg(mean[0]), round)*3600, label="$E(\omega_x^F)$");
# fig40_ax1.plot(time, np.round(np.rad2deg(flf[0]), round)*3600, label="$\omega_x^F$", alpha = 0.5, linestyle = "--", marker = "*");
fig40_ax1.grid(True)
fig40_ax1.legend(loc="lower center")
#дрейф северного акселерометра
fig40_ax2.set_title("Оценка скоростей\nдрейфа акселерометра\nканала Y")
fig40_ax2.set_xlabel("Мин")
# fig40_ax2.set_ylabel("Град/час")
# fig40_ax2.set_ylim(-0.25, 0.25)
fig40_ax2.plot(time, np.round(estimations["d_a_y"]), round, label="$\delta a_y$");
# fig40_ax2.plot(time, np.round(np.rad2deg(mean[1]), round)*3600, label="$E(\omega_y^F)$");
# fig40_ax2.plot(time, np.round(np.rad2deg(flf[1]), round)*3600, label="$\omega_y^F$", alpha = 0.5, linestyle = "--", marker = "*");
fig40_ax2.grid(True)
fig40_ax2.legend(loc="lower center",)
#дрейф вертикального акселерометра
# '''
fig40_ax3.set_title("Оценка скоростей\nдрейфа гироскопа\nканала Z")
fig40_ax3.set_xlabel("Мин")
# fig40_ax3.set_ylabel("Град/час")
# fig40_ax3.set_ylim(-0.25, 0.25)
fig40_ax3.plot(time, np.round(estimations["d_a_z"], round), label="$\delta a_z$");
fig40_ax3.grid(True)
fig40_ax3.legend(loc="lower center")
# '''
if args.SaveDir:
    fig40_path = os.path.join(args.SaveDir, "Оценки_дрейфов_акселерометров")
    saveFigure(fig40, fig40_path, format_type='pdf' if not args.PGF else "pgf")

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