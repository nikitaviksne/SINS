import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("/home/nikita/Документы/C_Cpp_progs/InertialNavigation/data/BLH_to_XYZ.csv", delimiter=";");

fig1, (fig1_ax1, fig1_ax2, fig1_ax3) = plt.subplots(3,1)
fig1_ax1.plot(data["X"], label="X")
fig1_ax1.legend(loc="best")
fig1_ax2.plot(data["Y"], label="Y")
fig1_ax2.legend(loc="best")
fig1_ax3.plot(data["Z"], label="Z")
fig1_ax3.legend(loc="best")

plt.show()
