import numpy as np

class R_Adaptive_Kalman(object):
    def __init__(self, dim_x, dim_z, x = None) -> None:
        self.dim_x = dim_x
        self.dim_z = dim_z
        # self.R = np.eye(dim_z)
        # self.Q = np.eye(dim_x)
        # self.P_aposteriori = np.eye(dim_x)
        self.P_aposteriori = np.zeros((dim_x, dim_x))
        self.x = x
        self.k_iter = 0
        self.C_k = 1e-5

    def predict(self):
        self.x_k_1 = self.Phi @ self.x # формирование априорной оценки
        self.P_apriori = self.Phi @ self.P_aposteriori @ self.Phi.T + (self.G @ self.Q) @ self.G.T
        pass
    
    def update(self, z):
        self.v_k = (z - self.H @ self.x_k_1).reshape(self.dim_z, 1) # расчет невязки
        self.C_k = self.k_iter/(self.k_iter + 1)*self.C_k + 1/(self.k_iter + 1)*(self.v_k @ self.v_k.T)
        self.R_k = self.C_k - self.H @ self.P_apriori @ self.H.T
        if (np.diag(self.R_k) < 0).any(): #проверяем на отрицаиельность диагональные элементы матрицы R
            # self.R_k = np.zeros((self.dim_z, self.dim_z))
            for i in range(self.dim_z):# и устанавливаем все диагональныеп элементы в нуль
                self.R_k[i,i] = 0
        # if np.linalg.det(self.H @ self.P_apriori @ self.H.T + self.R_k):
        #     self.k = (self.P_apriori @ self.H.T) @ np.linalg.inv(self.H @ self.P_apriori @ self.H.T + self.R_k)
        # else:
        self.k = (self.P_apriori @ self.H.T) @ np.linalg.pinv(self.H @ self.P_apriori @ self.H.T + self.R_k)
        
        self.P_aposteriori = (np.eye(self.dim_x) - self.k @ self.H) @ self.P_apriori
        self.x = self.x_k_1 + np.resize(self.k @ self.v_k, self.dim_x) # вычисление апостериорной оценки
        self.k_iter += 1
        pass
        # return self.x
