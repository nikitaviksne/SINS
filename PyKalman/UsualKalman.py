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
