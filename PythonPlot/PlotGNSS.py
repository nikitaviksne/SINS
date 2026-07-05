import argparse
import matplotlib.pyplot as plt
from pandas import read_csv as pd_read_csv

path = "/home/nikita/Документы/C_Cpp_progs/InertialNavigation/GNSS/Dubna/measurnments.csv"
data = pd_read_csv(path, delimiter=";");

PARAMETERS_PER_SAT = 10; #S,dS,X_i,Y_i,Z_i,tau_i, X_{i-1},Y_{i-1},Z_{i-1},tau_{i-1}

header = list(data)#.split(";")

fig1, (fig1_ax1, fig1_ax2) = plt.subplots(2,1)
#S
for s in [header[i] for i in range(1, len(header),PARAMETERS_PER_SAT) if not data[header[i]].empty]:
    fig1_ax1.plot(data[s], label=s)
#dS
for ds in [header[i] for i in range(2, len(header),PARAMETERS_PER_SAT)]:
    fig1_ax2.plot(data[ds], label=ds)
fig1_ax1.legend(loc='best')
fig1_ax2.legend(loc='best')

fig2, (fig2_ax1, fig2_ax2, fig2_ax3) = plt.subplots(3,1)
#X
for x in [header[i] for i in range(3, len(header),PARAMETERS_PER_SAT) if not data[header[i]].empty]:
    fig2_ax1.plot(data[x], label=x)
#Y
for y in [header[i] for i in range(4, len(header),PARAMETERS_PER_SAT)]:
    fig2_ax2.plot(data[y], label=y)
#Z
for y in [header[i] for i in range(5, len(header),PARAMETERS_PER_SAT)]:
    fig2_ax3.plot(data[y], label=y)
fig2_ax1.legend(loc='best')
fig2_ax2.legend(loc='best')
fig2_ax3.legend(loc='best')

plt.show()
