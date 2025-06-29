import numpy as np

class Kalman(object):
    def __init__(self, dim_x, dim_z, x = None) -> None:
        self.dim_x = dim_x
        self.dim_z = dim_z
        # self.R = np.eye(dim_z)
        # self.Q = np.eye(dim_x)
        # self.P_aposteriori = np.eye(dim_x)
        self.P_aposteriori = np.zeros((dim_x, dim_x))
        self.x = x

    def predict(self):
        self.P_apriori = self.Phi @ self.P_aposteriori @ self.Phi.T + (self.G @ self.Q) @ self.G.T
        self.k = (self.P_apriori @ self.H.T) @ np.linalg.inv(self.H @ self.P_apriori @ self.H.T + self.R)
        self.P_aposteriori = (np.eye(self.dim_x) - self.k @ self.H) @ self.P_apriori
        pass
    
    def update(self, z):
        # x_k = self.x
        x = np.dot(np.dot(self.H, self.Phi), self.x.reshape(self.dim_x, 1))
        self.x = (np.dot(self.Phi, self.x.T) + (self.k @ (z - x)).reshape(1, self.dim_x))[0]
        pass
        # return self.x

class Adaptive_Kalman(object):
    def __init__(self, dim_x, dim_z, x = None) -> None:
        self.dim_x = dim_x
        # self.dim_z = dim_z
        # self.R = np.eye(dim_z)
        # self.Q = np.eye(dim_x)
        # self.P_aposteriori = np.eye(dim_x)
        self.P_aposteriori = np.zeros((dim_x, dim_x))
        self.x = x
        self.k_iter = 0
        self.C_k = 1e-5

    def predict(self, z):
        self.x_k_1 = self.Phi @ self.x # формирование априорной оценки
        self.P_apriori = self.Phi @ self.P_aposteriori @ self.Phi.T + (self.G * self.Q) * self.G.T
        self.v_k = z - self.H @ self.x_k_1 # расчет невязки
        self.C_k = self.k_iter/(self.k_iter + 1)*self.C_k + 1/(self.k_iter + 1)*(self.v_k @ self.v_k.T)
        self.R_k = self.C_k - self.H @ self.P_apriori @ self.H.T
        if True in (self.R_k < 0):
            self.R_k = 0
        if np.linalg.det(self.H @ self.P_apriori @ self.H.T + self.R_k):
            self.k = (self.P_apriori @ self.H.T) @ np.linalg.inv(self.H @ self.P_apriori @ self.H.T + self.R_k)
        else:
            self.k = (self.P_apriori @ self.H.T) @ np.linalg.pinv(self.H @ self.P_apriori @ self.H.T + self.R_k)
        self.P_aposteriori = (np.eye(self.dim_x) - self.k @ self.H) @ self.P_apriori
        self.k_iter += 1
        pass
    
    def update(self):
        self.x = self.x_k_1 + self.k @ self.v_k # вычисление апостериорной оценки
        pass
        # return self.x
