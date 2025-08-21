import numpy as np
import quaternion

def norming(q:np.quaternion, delta_0:float, i:float):
    if np.abs(1-(q.w**2 + q.x**2 + q.y**2 + q.z**2)) > delta_0:
        q/=np.sqrt(np.abs(q.w**2 + q.x**2 + q.y**2 + q.z**2))
    return q

def matrix_to_quat(C:np.array):
    q_0 = 0.5*np.sqrt(C[0,0]+ C[1,1] + C[2,2] +1)
    if (q_0!=0):
        q_1 = (C[2,1] - C[1,2])/(4*q_0)
        q_2 = (C[0,2] - C[2,0])/(4*q_0)
        q_3 = (C[1,0] - C[0,1])/(4*q_0)
    else:
        q_1 = np.sqrt((C[0,0]+1)/2)
        q_2 = np.sqrt((C[1,1]+1)/2)
        q_3 = np.sqrt((C[2,2]+1)/2)
    q = np.quaternion(q_0, q_1, q_2, q_3)
    return q
def quat_to_matrix(Q:np.quaternion):
    c_11 = Q.w**2 + Q.x**2 - Q.y**2 - Q.z**2
    c_12 = 2*(Q.x*Q.y - Q.w*Q.z)
    c_13 = 2*(Q.x*Q.z + Q.w*Q.y)
    c_21 = 2*(Q.x*Q.y + Q.w*Q.z)
    c_22 = Q.w**2 + Q.y**2 - Q.x**2 - Q.z**2
    c_23 = 2*(Q.y*Q.z - Q.w*Q.x)
    c_31 = 2*(Q.x*Q.z - Q.w*Q.y)
    c_32 = 2*(Q.y*Q.z + Q.w*Q.x)
    c_33 = Q.w**2 + Q.z**2 - Q.x**2 - Q.y**2
    
    c_0 = np.sqrt(c_31**2 + c_33**2)

    C = np.array([[c_11, c_12, c_13], [c_21, c_22, c_23], [c_31, c_32, c_33]])
    return C

def quat_to_euler(Q: np.quaternion):
    c_11 = Q.w**2 + Q.x**2 - Q.y**2 - Q.z**2
    c_12 = 2*(Q.x*Q.y - Q.w*Q.z)
    c_13 = 2*(Q.x*Q.z + Q.w*Q.y)
    c_21 = 2*(Q.x*Q.y + Q.w*Q.z)
    c_22 = Q.w**2 + Q.y**2 - Q.x**2 - Q.z**2
    c_23 = 2*(Q.y*Q.z - Q.w*Q.x)
    c_31 = 2*(Q.x*Q.z - Q.w*Q.y)
    c_32 = 2*(Q.y*Q.z + Q.w*Q.x)
    c_33 = Q.w**2 + Q.z**2 - Q.x**2 - Q.y**2
        
    c_0 = np.sqrt(c_31**2 + c_33**2)

    C_b_ll = np.array([[c_11, c_12, c_13], [c_21, c_22, c_23], [c_31, c_32, c_33]])
        
    psi = np.arctan2(c_12, c_22)#*rad_to_deg
    gamma= - np.arctan2(c_31, c_33)#*rad_to_deg
    theta = np.arctan2(c_32, c_0)#*rad_to_deg

    C_b_ll = np.array([[c_11, c_12, c_13], [c_21, c_22, c_23], [c_31, c_32, c_33]])

    # psi = torch.arctan(c_12/c_22)       # in radians
    # gamma= - torch.arctan(c_31/c_33)    # in radians
    # theta = torch.arctan(c_32/c_0)      # in radians
    # psi = np.arctan(C_b_ll[0,1]/C_b_ll[1,1])                                     # in radians
    # gamma= - np.arctan(C_b_ll[2,0]/C_b_ll[2,2])                                  # in radians
    # theta = np.arctan(C_b_ll[2,1]/np.sqrt(C_b_ll[2,0]**2 + C_b_ll[2,2]**2))      # in radians

    angles = np.array([[psi],
                       [gamma], 
                       [theta]])

    return angles, C_b_ll

C_b_ll = lambda Psi, gamma, theta: np.array([[np.cos(gamma)*np.cos(Psi)+np.sin(gamma)*np.sin(theta)*np.sin(Psi), np.cos(theta)*np.sin(Psi), np.sin(gamma)*np.cos(Psi)-np.cos(gamma)*np.sin(Psi)*np.sin(theta)],
                    [-np.cos(gamma)*np.sin(Psi)+np.sin(gamma)*np.sin(theta)*np.cos(Psi), np.cos(theta)*np.cos(Psi), -(np.sin(gamma)*np.sin(Psi)+np.cos(gamma)*np.sin(theta)*np.cos(Psi))],
                    [-np.sin(gamma)*np.cos(theta), np.sin(theta), np.cos(gamma)*np.cos(theta)]]) #функция для матрицы ориентации из body в local level

g = 9.81
#проверка матрицы ориентации
psi = np.deg2rad(12);
gamma = np.deg2rad(32);
theta = np.deg2rad(2);

#C_o_b = np.array([0.0010410645, 1.0827141153, -0.0000124948, -1.0828304627, 0.0010397894, 0.1201860011, -0.0001411344, 0.0000114035, 0.9999082174]).reshape(3,3).T #все погрешности включены
# C_o_b = np.array([0.0000000000, 1.0839344336, 0.0000000000, -1.0840681197, 0.0000000000, 0.1205236714, -0.0001519857, 0.0000000000, 0.9998935784]).reshape(3,3).T #все погрешности выключены
C_o_b = C_b_ll(Psi=psi,gamma=gamma, theta=theta)
print(f"det(C) = {np.linalg.det(C_o_b)}")
print(f"matrix C_o_b is:\n{C_o_b}")
Ao = np.array([0,0,g]).reshape(3,1)
Ab = C_o_b @ Ao
print(f"Ab = {Ab}")
Q = matrix_to_quat(C_o_b);
print(f"Matrix to quaternion. Quaternion is\n{Q}")
C_o_b_q = quat_to_matrix(Q)
print(f"Quaternion to matrix. matrix is\n{C_o_b_q}")
Ab_q = C_o_b @ Ao
print(f"Ab_q = {Ab_q}")

angle, _ = quat_to_euler(Q);
print(np.rad2deg(angle))
exit(0)