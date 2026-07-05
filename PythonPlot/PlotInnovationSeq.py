'''
Считается невязки по псевдодальностям (измеренное минус расчетное) с учетом вращения Земли
'''
import matplotlib.pyplot as plt
# from pandas import read_csv as pd_read_csv
import pandas as pd
import numpy as np

U = np.deg2rad(15)/3600 # скорость вращения Земли рад/с
c = 299792458 # скорость света м/с 

def CalcR(CoordSat, Q):
    R1 = np.sqrt(pow(CoordSat[0] - Q[0], 2) + pow(CoordSat[1] - Q[1],2)+ pow(CoordSat[2] - Q[2],2));
    temp = U / c;
    R1 += 1*temp*(Q[1]*CoordSat[0] - Q[0]*CoordSat[1]);
    R1 += 1*Q[3]; #уход шкалы расчетное
    R1 += 1*CoordSat[3]; #уход шкалы спутника
    return R1;

path = "/home/nikita/Документы/C_Cpp_progs/InertialNavigation/GNSS/Dubna/measurnments_7Sat.csv"
data = pd.read_csv(path, delimiter=";");

PARAMETERS_PER_SAT = 10; #S,dS,X_i,Y_i,Z_i,tau_i, X_{i-1},Y_{i-1},Z_{i-1},tau_{i-1}

Q = [2797600.129, 2115824.843, 5309346.947, 133816.83808881836] # считаем, что это опорные координаты

header = list(data)#.split(";")


fig1, (fig1_ax1, fig1_ax2, fig1_ax3) = plt.subplots(3,1)
#S
for s in [header[i] for i in range(1, len(header),PARAMETERS_PER_SAT) if not data[header[i]].empty]:
    fig1_ax1.plot(data[s], label=s)
#dS
for ds in [header[i] for i in range(2, len(header),PARAMETERS_PER_SAT)]:
    fig1_ax2.plot(data[ds], label=ds)
# Z для 7-го спутника
Z = []# невязки
for i, row in data.iterrows():
    XYZSat = pd.Series([row["X_i(7)"], row["Y_i(7)"], row["Z_i(7)"], 1*row["Tau_i(7)"]])
    if not XYZSat.isnull().any():
        D = CalcR(XYZSat, Q)
        Z.append(row["S(7)"] - D)
        if (len(Z)>1):
            Q[3] -= Z[-1] - Z[-2]
fig1_ax3.plot(Z)
fig1_ax1.legend(loc='best')
fig1_ax2.legend(loc='best')

plt.show();


print(f"full = {len(data)}")
print(f"Z = {len(Z)}")