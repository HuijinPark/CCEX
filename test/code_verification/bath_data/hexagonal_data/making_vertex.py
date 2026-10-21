#!/usr/bin/env python 

import numpy as np

def Findvertex3D(CellShape,defect):
    Def = np.array(CellShape).T@np.array(defect)
    print("Def :",Def)

    vertex=[]
    for z in [0,1]:
        for y in [0,1]:
            for x in [0,1]:
                temp = np.array(CellShape).T@np.array([x,y,z])
                temp -= Def
                vertex.append(temp)
    
    return vertex


if __name__ == "__main__":
    
    CellShape=[[7.4619090824000001,  -12.9244056521000008,    0.0000000000000000],
     [7.4619090824000001,   12.9244056521000008,    0.0000000000000000],
     [0.0000000000000000,    0.0000000000000000,   12.9145771254999993]]

    defect = [0.3888888888991744, 0.444444444398421, 0.3750000000029061]

    vertex = Findvertex3D(CellShape, defect)
    
    print()
    for i in range(8):
        print("v%d :"%(i+1),vertex[i][0],vertex[i][1],vertex[i][2] )

