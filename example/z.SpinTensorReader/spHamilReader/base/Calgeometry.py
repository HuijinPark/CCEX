import numpy as np

#Calculation the distance between spins
def distance(I_pos=[],NV_info=[]):
    temp = 0.0
    for j in range(len(I_pos)):
        temp = temp + (I_pos[j]-NV_info[j])**2
    dis_I_S = np.sqrt(temp)
    return dis_I_S

#Picking out the spins within bath radius
def select(I_pos,NV_info,bath_r):
    P1bathwiRbath = []
    for i in range(len(I_pos)):
        if (distance(I_pos[i],NV_info) < bath_r):
            P1bathwiRbath.append(I_pos[i])
    return P1bathwiRbath

#Calculating the coefficient in hyperfine

def c_t(NV_info,align_I_position,dis_I_S):
    r_vec =align_I_position[2]- NV_info[2]
    return r_vec/dis_I_S
    
def s_t(NV_info,align_I_position,dis_I_S):
    r_vec = np.sqrt(((align_I_position[0]-(NV_info[0]))**2) + (((align_I_position[1]-NV_info[1])**2)))
    return r_vec/dis_I_S
    
def c_phi(NV_info,align_I_position):
    r_vec = np.sqrt(((NV_info[0]-align_I_position[0])**2) + ((NV_info[1]-align_I_position[1])**2))
    if r_vec == 0:
        return 1
    else:
        return (align_I_position[0]-NV_info[0])/r_vec
    
def s_phi(NV_info,align_I_position):
    r_vec = np.sqrt(((align_I_position[0]-NV_info[0])**2) + ((align_I_position[1]-NV_info[1])**2))
    if r_vec == 0:
        return 0
    else:
        return (align_I_position[1]-NV_info[1])/r_vec
