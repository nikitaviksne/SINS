import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("/home/nikita/Документы/C_Cpp_progs/InertialNavigation/data/debugKalman.csv", delimiter=";");

fig1, (fig1_ax1) = plt.subplots(1,1)
fig1_ax1.set_title("УМШВ")
#fig1_ax1.plot(data["Hx[umshv]"], label="Hx")
for i in range(1,8):
    fig1_ax1.plot(data[f"Z_{i}[umshv]"], label=f"Z_{i}")
    
fig1_ax1.legend(loc="best")

plt.show()
