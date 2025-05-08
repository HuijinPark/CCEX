import numpy as np

def print_atoms_from_unitcell(cell):
	for k in (cell.atoms).keys():
		for i in range(len(cell.atoms[k])):
			print(k,cell.atoms[k][i])

def print_atoms_from_BathArray(supercell,key=None):
	for i in range(len(supercell)):
		print(f"{i} )",supercell[i]['N'],supercell[i]['xyz'],end=" ")
		if key != None :
			if np.any(supercell[i][key]):
				print("!")
			print(supercell[i][key])
