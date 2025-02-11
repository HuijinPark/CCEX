from spHamilReader.base.rwfiles import readfile
import numpy as np

def checkLine(line,keyword):
	# function : 
	# 	check the line include all keywords
	boolList=[]
	for k in keyword.strip().split():
		boolList.append(k in line)
	return all(boolList)

def checkLine_all(line,keyword):
	# function : 
	# 	check the line is the same to keywords
	return (line == keyword.strip().split())

def checkTF(line,keyword):
	if 'T' in line[2]:
		return True
	else:
		return False

def readINCAR_spinHamil(fin):
	data = readfile(fin)
	opt = {}
	keys = ['LDMATRIX', 'LHYPERFINE', 'LEFG']
	for i in range(len(data)):
		try:
			if data[i][0] in keys:
				opt[data[i][0]] = checkTF(data[i],data[i][0])
		except:pass
	return opt

def readOUTCAR(fin,readincar=False):
	"""
		It extract data for last iteraction
		Args: 
			fin(list/str)   : data or file name
			readincar(bool) : Read or not INCAR part in OUTCAR 
		Return: 
			lastIterdata(list) : Last scf iteraction part ~ end of file
			incardata(dict) : incar part
			poscardata(dict) : poscar data - atomic species, total atomc #
	"""

	data = None
	try:
		data = readfile(fin)
	except:
		data = fin	

	# last iteration of outcar
	lastIterdata = []
	
	for i in range(len(data)-1,-1,-1):
		if data[i] != []:
			lastIterdata.append(data[i])
		if checkLine(data[i],"-"*39 + " Iteration "+"-"*39):
			break;

	# read poscar from outcar ( atomic species )
	poscardata = {"atom_data": [],
				  "unit_cell_information": [],
				  "total_atom_num": 0}	
	line = 0 
	while(line<len(data)):
		if checkLine(data[line], 'POTCAR:'):
			if data[line][2] not in poscardata["atom_data"]:
				poscardata["atom_data"].append(data[line][2])
			else:
				break;
		line += 1
	
	# read poscar from outcar ( the number of total atoms )
#	print("POSCAR : atomic position")
	while(line<len(data)):
		if checkLine_all(data[line], 'Primitive cell'):
			goread=False	
			while(not checkLine_all(data[line], 'ion indices of the primitive-cell ions')):

				if checkLine_all(data[line], 'position of ions in fractional coordinates (direct lattice)'):
					goread=True
					line+=1
				elif goread==True and data[line]==[]:
					goread=False

				if goread == True:
					poscardata["unit_cell_information"].append(list(map(float,data[line])))
					poscardata['total_atom_num'] +=1

				line+=1	
			break;
		line += 1

	# read incar from outcar 
	if readincar == True:
		incardata = {}	
		line = 0 
		while(line<len(data)):
			if checkLine(data[line], 'INCAR:'):
				line +=1
				while (data[line]!=[]):
					incardata[data[line][0]] = ' '.join(map(str,data[line][2:]))		
					line +=1
				break;
			line +=1

	if readincar == True:
		return (lastIterdata[::-1],poscardata, incardata)
	else:
		return (lastIterdata[::-1],poscardata)

if __name__ == '__main__':

	############################################
	#Read INCAR 
	#Display the options
	#opts = readINCAR_spinHamil("./src/strn90.in")
	#print(opts)

	############################################
	#Read OUTfile of QE data &
	#Extract the data for last iteration
	#In here, we will reversely read lines from the end of file
	outs,atoms,incar = readOUTCAR("./src/hfn100.out")

	from dictproc import print_dict
	#for j in range(len(outs)):
	#	print(outs[j])

	print_dict(incar)
	print_dict(atoms)
	############################################
