import re
import numpy as np
from spHamilReader.base.dictproc import *
import sys

#%% read POSCAR file
####################################################################
def read_POSCAR(inputFile,print_io=None):

    if type(inputFile) != str:
        sys.exit(f"{__name__} Error!")        

    if print_io!=None:
        print('\tRead POSCAR file - {}'.format(inputFile))
        print("\t","-"*50) 
    atomicData = {}

    latt_constant = []
    unit_cell_parameter = []
    atom_data = []
    atom_data_num = []
    unit_cell_information = []
    cell_tag=[]

    temp_unit_cell_information=[]
    with open(inputFile, 'r') as f:
        for i, line in enumerate(f):
            inputdata=line.strip().split()
            if i==1:
                latt_constant=list(map(float,inputdata))[0]
            elif i>=2 and i<=4:
                unit_cell_parameter.append(list(map(float,inputdata)))
            elif i==5:
                atom_data = inputdata
            elif i==6:
                atom_data_num = list(map(int,inputdata))
            elif i>=7:
                #atomic positions
                try:
                    unit_list=list(map(float,[inputdata[0],inputdata[1],inputdata[2]]))
                    if len(inputdata) == 3:
                        temp_unit_cell_information.append(unit_list)
                    elif len(inputdata) == 6:
                        unit_list_TF=unit_list+[inputdata[3],inputdata[4],inputdata[5]]
                        temp_unit_cell_information.append(unit_list_TF)
                    else:
                        sys.exit("Read_file error. check the selective dyanamics")    
                except:
                    if inputdata != []:
                        cell_tag.append(inputdata)

    f.close()

    num=0
    for i in range(len(atom_data)):
        if (i==0):num=0
        elif (i>0):num+=atom_data_num[i-1]
        for j in range(atom_data_num[i]):
            templist = []
            templist.append(atom_data[i])
            for k in range(len(temp_unit_cell_information[num+j])):
                templist.append(temp_unit_cell_information[num+j][k])
            unit_cell_information.append(templist)

    #atom_data=sum(atom_data,[])
#    tot_natm=sum(atom_data_num)

#    for j in range(0,len(atom_data)):
#        for i in range(0,atom_data_num[j]):
#            temp_list = []
#            temp_list.append(atom_data[j])
#            for k in range(len(temp_unit_cell_information[i])):
#                temp_list.append(temp_unit_cell_information[i][k])
#        unit_cell_information.append(temp_list)


    atomicData["latt_constant"] = latt_constant    
    atomicData["unit_cell_parameter"] = unit_cell_parameter
    atomicData["atom_data"] = atom_data
    atomicData["atom_data_num"] = atom_data_num
    atomicData["unit_cell_information"] = unit_cell_information
    atomicData["cell_tag"] = cell_tag 

    if print_io!=None:
        print_dict(atomicData)
        print("\t","-"*50) 
    return atomicData
####################################################################


#%% write new POSCAR file
####################################################################
def write_POSCAR(outputFile,new_atomicData,print_io=None): 
    # key : toprt  - To print information for the 
    #                    difference between old and new one,
    #                                    add old atomic data
    print('\n\tOutput POSCAR File : ',outputFile)
    print("\t","-"*50) 

    with open(outputFile, 'w') as f:
        f.write("POSCAR\n")
        f.write("{:>f}\n".format(new_atomicData["latt_constant"]))

        for i in range(len(new_atomicData["unit_cell_parameter"])):
            line="{:>15.10f}    {:>15.10f}    {:>15.10f}\n".format(
                float(new_atomicData["unit_cell_parameter"][i][0]),
                float(new_atomicData["unit_cell_parameter"][i][1]),
                float(new_atomicData["unit_cell_parameter"][i][2]))
            f.write(line)
        
        for i in range(len(new_atomicData["atom_data"])):
            f.write("{} ".format(new_atomicData["atom_data"][i]))
        f.write("\n")

        for i in range(len(new_atomicData["atom_data_num"])):
            f.write("{} ".format(new_atomicData["atom_data_num"][i]))
        f.write("\n")

        for i in range(len(new_atomicData["cell_tag"])):
            for j in range(len(new_atomicData["cell_tag"][i])):        
                f.write("{} ".format(new_atomicData["cell_tag"][i][j]))
            f.write("\n")
                
        for i in range(len(new_atomicData["unit_cell_information"])):
            for j in range(1,len((new_atomicData["unit_cell_information"][i]))):
                if j < 4: 
                    line="{:>15.10f}   ".format(new_atomicData["unit_cell_information"][i][j])
                    f.write(line)
                else:
                    try:
                        line="{:>d} ".format(int(new_atomicData["unit_cell_information"][i][j]))
                        f.write(line)
                    except:
                        line="{:>s} ".format(new_atomicData["unit_cell_information"][i][j])
                        f.write(line)
            f.write("\n")
    f.close()       

    if print_io!=None:
        print_dict(new_atomicData)
        print("\t","-"*50) 

    print("\t","-"*50,"\n") 
####################################################################

if __name__ =="__main__":
    import dictproc as dp

    f = "/home/huijin/cal/vasp/10_hBN_gStrain{hsym}/results_20A/contcars/POSCAR_g90_VB_S1"
    atomicData = read_POSCAR(f)
    for i,v in enumerate(atomicData["unit_cell_information"][-10:]):
        print(i,v)

