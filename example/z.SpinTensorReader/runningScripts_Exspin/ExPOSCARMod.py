from ioPoscar import read_POSCAR, write_POSCAR
import sys


fin = sys.argv[1]
fout = sys.argv[2]
Exspin = sys.argv[3]

data=read_POSCAR(fin)

if Exspin=='NVHn_up':

    data["atom_data"] += ['VC']
    data["atom_data"][1] = 'H'
    data["atom_data"][3] = 'C'
    data["atom_data_num"] += [1]
    data["unit_cell_information"]+=[['VC',0.5,0.5,0.5]]

elif Exspin=='NVHn_dw':

    data["atom_data"] += ['VC']
    data["atom_data"][1] = 'H'
    data["atom_data"][3] = 'C'
    data["atom_data_num"] += [1]
    data["unit_cell_information"]+=[['VC',0.5,0.5,0.625]]

if Exspin=='NVHn_up2':

    data["atom_data"] += ['VC']
    data["atom_data"][1] = 'H'
    data["atom_data"][2] = 'C'
    data["atom_data_num"] += [1]
    data["unit_cell_information"]+=[['VC',0.5,0.5,0.5]]

elif Exspin=='NVHn_dw2':

    data["atom_data"] += ['VC']
    data["atom_data"][1] = 'H'
    data["atom_data"][2] = 'C'
    data["atom_data_num"] += [1]
    data["unit_cell_information"]+=[['VC',0.5,0.5,0.625]]

elif Exspin=='N2p_up':

    data["atom_data"] = ['C', 'N', 'C', 'ele']
    data["atom_data_num"] = [568, 2, 6, 1 ]
    data["unit_cell_information"]+=[data["unit_cell_information"][-8]]

elif Exspin=='NV0_up':

    data["atom_data"] = ['C', 'N', 'C', 'VC']
    data["atom_data_num"] = [571, 1, 3, 1 ]
    data["unit_cell_information"]+=[['VC',0.5,0.5,0.5]]

elif Exspin=='NV0_dw':

    data["atom_data"] = ['C', 'N', 'C', 'VC']
    data["atom_data_num"] = [571, 1, 3, 1 ]
    data["unit_cell_information"]+=[['VC',0.5,0.5,0.625]]

elif Exspin=='NVn_up':

    data["atom_data"] = ['C', 'N', 'C', 'VC']
    data["atom_data_num"] = [571, 1, 3, 1 ]
    data["unit_cell_information"]+=[['VC',0.5,0.5,0.5]]

elif Exspin=='NVn_dw':

    data["atom_data"] = ['C', 'N', 'C', 'VC']
    data["atom_data_num"] = [571, 1, 3, 1 ]
    data["unit_cell_information"]+=[['VC',0.5,0.5,0.625]]


elif Exspin=='VH0_up':

    data["atom_data"] = ['C', 'H', 'C', 'VC']
    data["atom_data_num"] = [572, 1, 3, 1 ]
    data["unit_cell_information"]+=[['VC',0.5,0.5,0.5]]

elif Exspin=='VH0_dw':

    data["atom_data"] = ['C', 'H', 'C', 'VC']
    data["atom_data_num"] = [572, 1, 3, 1 ]
    data["unit_cell_information"]+=[['VC',0.5,0.5,0.625]]

elif Exspin=='VHn_up':

    data["atom_data"] = ['C', 'H', 'C', 'VC']
    data["atom_data_num"] = [572, 1, 3, 1 ]
    data["unit_cell_information"]+=[['VC',0.5,0.5,0.5]]

elif Exspin=='VHn_dw':

    data["atom_data"] = ['C', 'H', 'C', 'VC']
    data["atom_data_num"] = [572, 1, 3, 1 ]
    data["unit_cell_information"]+=[['VC',0.5,0.5,0.625]]

elif Exspin=='Vn_up':

    data["atom_data"] = ['C', 'C', 'VC']
    data["atom_data_num"] = [571, 4, 1 ]
    data["unit_cell_information"]+=[['VC',0.5,0.5,0.5]]

elif Exspin=='Vn_dw':

    data["atom_data"] = ['C', 'C', 'VC']
    data["atom_data_num"] = [571, 4, 1 ]
    data["unit_cell_information"]+=[['VC',0.5,0.5,0.625]]

elif Exspin=='P1_up':

    data["atom_data"] = ['C', 'N', 'C', 'ele']
    data["atom_data_num"] = [571, 1, 4,1 ]
    data["unit_cell_information"]+=[data["unit_cell_information"][-5]]

else:
    sys.exit(f"There are no options : {Exspin}")
write_POSCAR(fout,data)

