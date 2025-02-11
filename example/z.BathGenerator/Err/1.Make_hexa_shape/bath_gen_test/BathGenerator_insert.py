#!/usr/bin/env python
#-*- coding:utf-8 -*-

import random as r
import numpy as np
import numpy.linalg as lin
import os
import sys
import math
import re
import time
from collections import Counter

#How to use this code?
###############################################################
#
#1. module load anaconda3
#
#2-1. ./BathGenerator_insert.py [periodic inputfile] [Defect inputfile] [probability file] [size] [Impurity atom Name] [Vacancy atom Name] [outfile name] [inital number] [range number]
#
#--> only make the [size, size, size] cubic outfile using arbitray structure(POSCAR)
#--> if you don't want to make vaccancy atom, write the Vacancy atom Name = "NONE"
#
#
#2-2. ./BathGenerator_insert.py [periodic inputfile] [Defect inputfile] [probability file] [size] [Impurity atom Name] [Vacancy atom Name] [outfile name] [inital number] [range number] [2D option]
#--> [2D option] = '2D' or '2' (including "2")
#
#
###############################################################
#periodic inputfile     :   POSCAR which has bulk or pure structure
#Defect inputfile       :   POSCAR or csv file having the defect atom
#outputfile             :   the name you wanted
#size                   :   the outfile's bath size if you want to [50A,50A,50A],
#                           size is 50
#density                :   [ppm]=[10^-4%] --> 1.1% = 11000ppm
#inital number          :   write the "1"
#range number           :   you want to get the only one file, 
#                           the range number is 1.
#                           or you want a 3 files, write the 3
###############################################################



#1. Function define
#read the inputfile data and find the unit cell information
############################
def unitCellInfo(inputfile):
    print('File name is : ',inputfile)
    with open(inputfile, 'r') as f:
        
        unit_cell_parameter=[]
        atom_data=[]
        atom_data_num=[]
        unit_cell_information=[]
        temp_unit_cell_information=[]

        for i, line in enumerate(f):
            inputdata=line.strip().split()
            if i == 1:
                LatMult=list(map(float,inputdata))[0]
            elif i >= 2  and i <= 4:
                unit_cell_parameter.append(list(map(float,inputdata)))
            elif i == 5:
                atom_data.append(inputdata)
            elif i == 6:
                atom_data_num.append(list(map(int,inputdata)))
            elif i > 7:
                unit_list=list(map(float,inputdata))
#                for j in range(len(unit_list)):
#                    test_one=unit_list[j]
#                    if test_one >= 1:
#                        unit_list[j] = round(test_one-1,10)
                temp_unit_cell_information.append(unit_list)
    f.close()

    #multitple unit_cell_parameter
    #print("LatMult : ",LatMult)
    for i in range(len(unit_cell_parameter)):
        for j in range(3):
            unit_cell_parameter[i][j]=unit_cell_parameter[i][j]*LatMult

    #form 2D array to 1D array
    atom_data=sum(atom_data,[])
    atom_data_num=sum(atom_data_num,[])

    if len(atom_data) != len(atom_data_num):
        print("plz, check the atom species and number!!\n")
        sys.exit()
    elif len(unit_cell_parameter) != 3:
        print("plz, check the cell parameter!!\n")
        sys.exit()
    elif len(temp_unit_cell_information) != sum(atom_data_num):
        print("plz, check the cell informatioin number!!\n")
        sys.exit()

    for j in range(0, len(atom_data)):
        for i in range(0, atom_data_num[j]):
            temp_list=[atom_data[j], temp_unit_cell_information[0][0],
                        temp_unit_cell_information[0][1],
                        temp_unit_cell_information[0][2]]
            del temp_unit_cell_information[0]
            unit_cell_information.append(temp_list)
    
    return unit_cell_parameter, atom_data, atom_data_num, unit_cell_information;
        

############################
def dist(pos1,pos2):    #array type
    tmp = (pos2 - pos1)
    return np.math.sqrt((tmp[0])**2 + (tmp[1])**2 + (tmp[2])**2)


#2. Defect information -> remove the overlap one.
############################
def Defect_position(unit_cell_information, Impurity, Vacancy, latticeParameter,sizeX,sizeY,sizeZ,Opt_2D):

    #now the unit_cell_information is str type due to atomic name
    #And the atomic name is 
    atomic_name=np.array(unit_cell_information)
    atomic_name=atomic_name[:,0]


    #find the mid lattice numbering 
    #Mid_point=MatSolv(latticeParameter,size/2,size/2,size/2):
    Lattic=np.array(latticeParameter)
    LatticT=Lattic.T
    inv_Lat=lin.inv(LatticT)
    Position=[[sizeX/2],[sizeY/2],[sizeZ/2]] #angstrom
    Mid_point=np.dot(inv_Lat, Position)

    MidA=int(Mid_point[0])
    MidB=int(Mid_point[1])
    MidC=int(Mid_point[2])

    a_list = range(int(MidA-2)+1,math.ceil(MidA+2))
    b_list = range(int(MidB-2)+1,math.ceil(MidB+2))
    c_list = range(int(MidC-2)+1,math.ceil(MidC+2))

    print("a_list(in Defect) : ",a_list)
    print("b_list(in Defect) : ",b_list)
    print("c_list(in Defect) : ",c_list)

    ref_length = lin.norm(np.array(latticeParameter[2]))#a3 length (angstrom)

    #Defect_p    = np.array([0.5,0.5,0.5]) #set ref position
    Defect_p    =np.array([sizeX/2,sizeY/2,sizeZ/2]) #ref position

    #tmpDistance = ref_length # angstrom, if you can find defect --> increase
    #Caution #Caution #Caution
    tmpDistance = ref_length*3 # angstrom, if you can find defect --> increase
    #Caution #Caution #Caution
    n = np.array([0.0,0.0,0.0]) # to find close spin central
    Plus_lat=[0,0,ref_length/2] # add ref vacancy position (h_BN)

    #if the matrial is 2D, we need to check for the defect position in the 0-th layer
    p=re.compile("2")
    if (Opt_2D):
        #ref position
        Defect_p  =np.array([sizeX/2,sizeY/2,sizeZ/2]) 
        c_list=range(1)
        ref_length=lin.norm(np.array(latticeParameter[1]))
        Plus_lat=[0,ref_length/2,0] # add ref vacancy position

#    print("NV position candidate list :\n")
    for a in a_list:
        for b in b_list:
            for c in c_list:
                for p0 in range(0,len(unit_cell_information)):
                    if atomic_name[p0] == Impurity:
                        Data_tmp=[a,b,c,p0]
                        tmpN=RealPosition(latticeParameter,unit_cell_information,Data_tmp)

                        distance = dist(tmpN,Defect_p)
                        if distance <= tmpDistance:
                            tmpDistance = distance 
                            #print("tempDistance : ",tmpDistance)
                            N_unit_data=Data_tmp
                            n = tmpN
                            #print('test for tmpN',tmpN)
                            #print('unit_cell data ', unit_cell_information[p0])
    print("N real site      : ",n)
    print("N_unit_data      : ",N_unit_data )

    # V => spin position + zVector
    # |
    # N => spin position and this position should be existed in nearest 0.5, 0.5, 0.5
    # nvArr = [N, V]

    #v_tmp   = n + vectorZ/2 #for BZO, NV diamond
    if(Vacancy != 'NONE'):
        v_tmp   = n + Plus_lat #for z-axis (angstrom)
        checkDistance = ref_length*3 #angstrom
    
        print('origin vacancy       : ', v_tmp)

    #the vacancy is in the same lattice
        a_list=[N_unit_data[0]]
        b_list=[N_unit_data[1]]
        c_list=[N_unit_data[2]]
    # test the v site is real
        for a in a_list:
            for b in b_list:
                for c in c_list:
                    for p0 in range(0,len(unit_cell_information)):
                        if (atomic_name[p0] == Vacancy and p0 != N_unit_data[3]) :
                            Data_tmp1=[a,b,c,p0]
                            tmpV=RealPosition(latticeParameter,unit_cell_information,Data_tmp1)

                            distance = dist(tmpV, v_tmp)
                            if (distance <= checkDistance) :
                                checkDistance = distance
                                v = tmpV
                                V_unit_data=Data_tmp1
        print("V real site      : ",v)
        print("V_unit_data      : ",V_unit_data )


        return [n,v], [N_unit_data, V_unit_data]

    else:
        return [n], [N_unit_data]

#3. read the nucelar spin probability of isotopes nuclear
############################
def probability(concTempfile, atom_data):
    with open(concTempfile,'r') as f:
        pro_name=[] #Atom name
        pro_spin=[] #Nuclear spin of atom
        pro_num=[]  #probability of nuclear spin and [%] unit

        for m, line in enumerate(f):
            data=line.split()
            if m % 3 == 0:
                pro_name.append(data)
            elif m % 3 == 1:
                pro_spin.append(data)
            elif m % 3 == 2:
                pro_num.append(list(map(float,data)))
        if(len(pro_num)+len(pro_spin)+len(pro_name))%3 != 0:
            print("plz, check the probability fiel and number error !!\n")
   
    #form 2D array to 1D array
    pro_name = sum(pro_name,[])
    pro_spin = sum(pro_spin,[])
    pro_num  = sum(pro_num,[])
    print("pro_name :",pro_name)
    print("pro_spin :",pro_spin)
    print("pro_num :",pro_num)
    

    f.close()
        
    return pro_name, pro_spin, pro_num

#4. solve the matrix to know the lattice index
############################
def MatSolv(latticeParameter,A,B,C):
#
#(i*a1 + j*b1 + k*c1, i*a2 + j*b2 + k*c2, i*a3 + j*b3 + k*c3)
# = (cubic X, cubic Y, cubic Z) = (A, B, C)

    #               a1, a2, a3      
    # [i j k]  dot  b1, b2, b3   =  [A B C]
    #               c1, c2, c3      
    #
    # [i j k] is # of a,b,c
    # [A B C] is length (unit : angstrom)
    #
    #--> knowing [i,j,k] ==> result
    Lattic=np.array(latticeParameter)
    LatticT=Lattic.T
    inv_Lat=lin.inv(LatticT)
    Position=[[A],[B],[C]]
    result=np.dot(inv_Lat, Position)

    return result

#5.real position of atom 
############################
def RealPosition(latticeParameter,unit_cell_information,data_set):
    LengthX     = np.array(latticeParameter[0])
    LengthY     = np.array(latticeParameter[1])
    LengthZ     = np.array(latticeParameter[2])

    A=data_set[0]
    B=data_set[1]
    C=data_set[2]
    p0=data_set[3]

    #(real site)
    Real = (unit_cell_information[p0][1]+A)*LengthX+(unit_cell_information[p0][2]+B)*LengthY+(unit_cell_information[p0][3]+C)*LengthZ      
    #print("unit_cell_information : ",unit_cell_information)
    return Real


#6.judge the boundary of cubic size
############################
def BDcut(Ispin_info_temp,sizeX,sizeY,sizeZ):
    
    maxlat=[sizeX,sizeY,sizeZ]
    minlat=[0,0,0] 

    ConMax=np.all(Ispin_info_temp <= maxlat)
    ConMin=np.all(Ispin_info_temp >= minlat)
    
    if ConMax and ConMin:
        return True
    else:
        return False

#7.judge the 2D option
############################
def ReadOption(Option2D):
    p=re.compile("2")
    if p.search(Option2D):
        return True
    else:
        return False

######################################################
######################################################
#main code start

#1. Input the information ############################
print('\n----------- System Option!! ----------\n')
print("system option length : ",len(sys.argv))
print("system option : \n",sys.argv)

if (len(sys.argv) != 10 and len(sys.argv) != 11 ):
    print("your option is poor!!")
    print("[1] : importFilename, [2] : insert file name, [3] : probability File, [4] : bath size(cubic), [5] : Impurity atom name, [6] : Vacancy atom name, [7] : save file name, [8] : inital file number, [9] : last file number, [10] : 2D option ")

inputfile       = sys.argv[1]
insertfile      = sys.argv[2]
concTempfile    = sys.argv[3]
size            = int(sys.argv[4])  #angstrom
                                    #make the lattice as cubic.
Impurity        = sys.argv[5]
Vacancy         = sys.argv[6]
outputfile      = sys.argv[7]
numOfConfigIni  = int(sys.argv[8])  #inital numbering
numOfConfig     = int(sys.argv[9])  #range output data file

if (len(sys.argv) == 11):
    print(len(sys.argv))
    Option2D       = sys.argv[10]
    Opt_2D=ReadOption(Option2D)
    print("2D_axis : ",Opt_2D)
else :
    Opt_2D=False
print("2D_axis option : ",Opt_2D)
    

#2. Other parameter ##################################
print('\n----------- Inputfile and Setting number!! ----------\n')
con_unit            = 0.01 #due to [unit : %] (1% = 1000000ppm)
                           #convert [% unit] to numbering
#about bulk file
print('Read bulk file!\n')
latticeParameter, atom_data, atom_data_num, unit_cell_information = unitCellInfo(inputfile)

#Insert file 
print('Read insert file!\n')
insert_latticeParameter, insert_atom_data, insert_atom_data_num, insert_unit_cell_information = unitCellInfo(insertfile)
print("insert_unit_cell_information : ", insert_unit_cell_information)

if (latticeParameter != insert_latticeParameter):
    print("Nope, the lattice parameter is not the same!!")
    sys.exit()


pro_name, pro_spin, pro_num = probability(concTempfile, atom_data)



#unitCellnumber : change the min and max using MatSolv
#cubic boundary : [0,size] 
XlatNumber=[0,0]
YlatNumber=[0,0]
ZlatNumber=[0,0]
sizeX=size #cubic X size
sizeY=size #cubic Y size
sizeZ=size #cubic Z size 

#if you want to make the 2D (no z-axis),
#to make the below for loop
if (Opt_2D):
    sizeZ=0

for i in 0, sizeX:
    for j in 0, sizeY:
        for k in 0, sizeZ:
            count_lat=MatSolv(latticeParameter,i,j,k)
            if count_lat[0] < XlatNumber[0]:  #max Xlat
                XlatNumber[0]=count_lat[0]
            elif count_lat[0] > XlatNumber[1]: #min Xlat
                XlatNumber[1]=count_lat[0]

            if count_lat[1] < YlatNumber[0]: #max Ylat
                YlatNumber[0]=count_lat[1]
            elif count_lat[1] > YlatNumber[1]: #min Ylat
                YlatNumber[1]=count_lat[1]

            if count_lat[2] < ZlatNumber[0]: #max Zlat
                ZlatNumber[0]=count_lat[2]
            elif count_lat[2] > ZlatNumber[1]: #min Zlat
                ZlatNumber[1]=count_lat[2]

unitCellnumber = [[int(XlatNumber[0]),int(XlatNumber[1])],
                [int(YlatNumber[0]),int(YlatNumber[1])],
                [int(ZlatNumber[0]),int(ZlatNumber[1])]]

#to make the one size lattice 
if (Opt_2D):
    sizeZ=latticeParameter[0][2]+latticeParameter[1][2]+latticeParameter[2][2]

print('Size of SuperCell (angstrom) : ',[sizeX,sizeY,sizeZ])
print('UnitCell # : ',unitCellnumber)

print('\nunitcell info : \n',unit_cell_information)
print('\nlatticeParameter       : \n',latticeParameter)

#Find the Defect center spin position
print('\n----------- Defect site finding!! ----------\n')
#Defect_info_tmp, Defect_data = Defect_position(unit_cell_information, Impurity, Vacancy,latticeParameter,sizeX,sizeY,sizeZ,Opt_2D)
print("We find the Defect position in the",insertfile)
Defect_info_tmp, Defect_data = Defect_position(insert_unit_cell_information, Impurity, Vacancy,latticeParameter,sizeX,sizeY,sizeZ,Opt_2D)
print('Length of Defect_information : ', len(Defect_info_tmp))

N_info_tmp = Defect_info_tmp[0]
N_data_tmp = Defect_data[0]
print('\nN real site        : ',N_info_tmp)
print('N data             : ',N_data_tmp)

#Impurity site site (real site)
Defect_info_N=RealPosition(latticeParameter,insert_unit_cell_information,N_data_tmp) 

print('Defect of impurity : ',Defect_info_N)
Defect_site_tmp=[]
Defect_site_tmp.append(N_info_tmp.tolist())


if (len(Defect_info_tmp) == 2):
    V_info_tmp = Defect_info_tmp[1]
    V_data_tmp = Defect_data[1]
    print('\nV real site        : ',V_info_tmp)
    print('V data             : ',V_data_tmp)

    #Vacancy site (real site)
    Defect_info_V=RealPosition(latticeParameter,insert_unit_cell_information,V_data_tmp)  

    print('Defect of vacancy  : ',Defect_info_V)
    Defect_site_tmp.append(V_info_tmp.tolist())

#show the all Defect site 
print('\nDefect_site_tmp : \n',Defect_site_tmp)


#setting the isotope atom number in unit cell
print('\n----------- Setting Isotope atom number!! ----------\n')
isotope_data_num=[]
isotope_data_name=[]
p=re.compile('[a-zA-Z]+')
for i in range(0,len(pro_num)):
    iso_name=p.findall(pro_name[i])
    if len(iso_name) != 1:
        print("The error of isotope name in probability file !!\n")
        print("Change the correct isotope atom name\n")
        sys.exit()
    elif(iso_name[0] not in atom_data) and (iso_name[0] not in insert_atom_data):
        if(iso_name[0] is not Impurity) or (iso_name[0] is not Vacancy):
            print("The atom name of isotope is not in the unit cell\n")
            print("plz, check the isotope file(probablity file)\n")
            print("Checking iso name :",iso_name[0])
            sys.exit()
        else:
            print("Is this the Vacancy or Impurity?\n")
            print("Checking iso name :", iso_name[0])
    else:
        if iso_name[0] in atom_data:
            index_num=atom_data.index(iso_name[0])
            isotope_data_num.append(atom_data_num[index_num])
            isotope_data_name.append(atom_data[index_num])
        elif iso_name[0] in insert_atom_data:
            index_num=insert_atom_data.index(iso_name[0])
            isotope_data_num.append(insert_atom_data_num[index_num])
            isotope_data_name.append(insert_atom_data[index_num])

print('Isotope data Name : ', isotope_data_name)
print('Isotope data Number : ', isotope_data_num)


#3. Random generator spin within supercell size
print('\n----------- Random generator!! ----------\n')

#atomic name
atomic_name=np.array(unit_cell_information)
atomic_name=atomic_name[:,0]

insert_atomic_name=np.array(insert_unit_cell_information)
insert_atomic_name=insert_atomic_name[:,0]

loopX=range(unitCellnumber[0][0], unitCellnumber[0][1]+1)
loopY=range(unitCellnumber[1][0], unitCellnumber[1][1]+1)
loopZ=range(unitCellnumber[2][0], unitCellnumber[2][1]+1)

if (Opt_2D):
    loopZ=range(1)

#checking for configure atom
Checking_atom=0

#print or not print detail data
#Onprint=True
Onprint=False
if Onprint==False:print("...\nSkip atom detail print.\nif you want to see detail, Change the option : 'Onprint=True'\n...\n")

#main code for generator of spin
for q in range(numOfConfigIni, numOfConfig+1):

    Ispin_info = []

    for i in loopX:
        for j in loopY:
            for k in loopZ:
                if Onprint:print("SuperCell[%d,%d,%d], Ispin_info_temp \n:"%(i,j,k))
                #About the defect atom unitcell
                if (i == N_data_tmp[0]) and (j==N_data_tmp[1]) and (k==N_data_tmp[2]):
                    for n in range(len(insert_unit_cell_information)):
                        Ispin_info_temp = []
                        Ispin_data_temp=[i,j,k,n]
                        Ispin_info_temp=RealPosition(latticeParameter,insert_unit_cell_information,Ispin_data_temp)
                        if Onprint:{print("Ispin data info : ",Ispin_data_temp),
                                    print("atomic name     : ",insert_atomic_name[n]),
                                    print("Ispin info      : ", Ispin_info_temp)}
                        if (Ispin_data_temp in Defect_data):
                            if Onprint:print("X(overlapped)\n")
                        else:
                            #example of BZO, atom_data : [Ba, Zr, O]
                            tot_pro_num=0
                            Iso_count=[]
                            p=re.compile(insert_atomic_name[n])
    
                            for iso in range(len(pro_num)):
                                if(p.search(pro_name[iso])!= None):
                                    tot_pro_num=tot_pro_num+pro_num[iso]
                                    Iso_count.append(iso)
                            if Onprint:print("Iso_count       : ",Iso_count)
    
                            #probability is between 0 and 1
                            prob = tot_pro_num *con_unit
                            if Onprint:print("prob            : ",prob)
    
                            #ref random value is between 0 and 1
                            randVal = r.uniform(0,1)
                            if Onprint:print("randVal         : ",randVal)
    
                            if (randVal <= prob):
                                #consider the Isotope atom is one
                                if (len(Iso_count) == 1):
                                    if Onprint:{print("n           :",n),
                                                print("pro_name[n] : ",pro_name[Iso_count[0]])}
                                    Ispin_info.append(np.hstack([Ispin_info_temp,pro_name[Iso_count[0]]]))
    
                                #consider the Isotope atom is more one
                                else:
                                    length_count = len(Iso_count)
                                    IsoP=[]
                                    for IC in range(length_count):
                                        IsoP.append(pro_num[Iso_count[IC]]/tot_pro_num)
                                    SetIso=np.random.choice(length_count,1,p=IsoP)
                                    #np.random.choice(list,# items, probability)
                                    #-->pick the item # in the list with probability
                                    if Onprint:{print("IsoP        : ",IsoP),
                                                print("n           : ",n),
                                                print("SetIso      : ",SetIso[0]),
                                                print("Iso_count[SetIso] : ",Iso_count[SetIso[0]]),
                                                print("pro_name[n] :",pro_name[Iso_count[SetIso[0]]])}
                                    Ispin_info.append(np.hstack([Ispin_info_temp,pro_name[Iso_count[SetIso[0]]]]))
    
                                if Onprint:print("O(rand)\n")
                                Checking_atom=Checking_atom+1
    
                            else:
                                if Onprint:print("X(rand)\n")
                                Checking_atom=Checking_atom+1



                #about the other unit cell atoms
                else:
                    for n in range(len(unit_cell_information)):
                        Ispin_info_temp = []
                        Ispin_data_temp=[i,j,k,n]
                        Ispin_info_temp=RealPosition(latticeParameter,unit_cell_information,Ispin_data_temp)
                        if Onprint:{print("Ispin data info : ",Ispin_data_temp),
                                    print("atomic name     : ",atomic_name[n]),
                                    print("Ispin info      : ", Ispin_info_temp)}
    
                        JudgeBD=BDcut(Ispin_info_temp,sizeX,sizeY,sizeZ)
                        if (Ispin_data_temp in Defect_data):
                            if Onprint:print("X(overlapped)\n")
                        elif not JudgeBD :
                            if Onprint:print("X(Boundary out)\n")
                        else:
                            #example of BZO, atom_data : [Ba, Zr, O]
                            tot_pro_num=0
                            Iso_count=[]
                            p=re.compile(atomic_name[n])
    
                            for iso in range(len(pro_num)):
                                if(p.search(pro_name[iso])!= None):
                                    tot_pro_num=tot_pro_num+pro_num[iso]
                                    Iso_count.append(iso)
                            if Onprint:print("Iso_count       : ",Iso_count)
    
                            #probability is between 0 and 1
                            prob = tot_pro_num *con_unit
                            if Onprint:print("prob            : ",prob)
    
                            #ref random value is between 0 and 1
                            randVal = r.uniform(0,1)
                            if Onprint:print("randVal         : ",randVal)
    
                            if (randVal <= prob):
                                #consider the Isotope atom is one
                                if (len(Iso_count) == 1):
                                    if Onprint:{print("n           :",n),
                                                print("pro_name[n] : ",pro_name[Iso_count[0]])}
                                    Ispin_info.append(np.hstack([Ispin_info_temp,pro_name[Iso_count[0]]]))
    
                                #consider the Isotope atom is more one
                                else:
                                    length_count = len(Iso_count)
                                    IsoP=[]
                                    for IC in range(length_count):
                                        IsoP.append(pro_num[Iso_count[IC]]/tot_pro_num)
                                    SetIso=np.random.choice(length_count,1,p=IsoP)
                                    #np.random.choice(list,# items, probability)
                                    #-->pick the item # in the list with probability
                                    if Onprint:{print("IsoP        : ",IsoP),
                                                print("n           : ",n),
                                                print("SetIso      : ",SetIso[0]),
                                                print("Iso_count[SetIso] : ",Iso_count[SetIso[0]]),
                                                print("pro_name[n] :",pro_name[Iso_count[SetIso[0]]])}
                                    Ispin_info.append(np.hstack([Ispin_info_temp,pro_name[Iso_count[SetIso[0]]]]))
    
                                if Onprint:print("O(rand)\n")
                                Checking_atom=Checking_atom+1
    
                            else:
                                if Onprint:print("X(rand)\n")
                                Checking_atom=Checking_atom+1
                if Onprint:print("\n")

#to prepare the 4-3
only_spin2=np.array(Ispin_info)

#to check the atom probablitiy 
print("---------------------------------------------")
print("<Ispin length infomation>")
Num_spin=Counter((np.array(Ispin_info)[:,3]))
print("Total atom #")
Num_allatom=0
for key in Num_spin:
    print(key,":",Num_spin[key])
    Num_allatom=Num_allatom+Num_spin[key]
print("\nAtom percentage[%] in only spin:")
for key in Num_spin:
    print(key,":",round(Num_spin[key]/Num_allatom*100,5),"%")
print("--> All spin number[#] :", Num_allatom)

print("\nAtom percentage[%] in configure (except for Defect atom)")
for key in Num_spin:
    print(key,":",round(Num_spin[key]/Checking_atom*100,6),"%")
print("--> All atom number[#] :",Checking_atom)
print("---------------------------------------------")


#4-1. Save the information ####################################
print("\nNow Save the information!-----------------------")
f = open(outputfile+"_%d"%(q), 'w')
only_spin=np.array(Ispin_info)
print(only_spin)

if(len(only_spin) != 0):

    only_spin=only_spin[:,3]
    Ispin_info=np.array(Ispin_info)
    Ispin_info=Ispin_info[:,:3].astype(np.float)

    line = "{:>10.5f}\t{:>10.5f}\t{:>10.5f}\t{:>10.5f}\n".format(float(len(Ispin_info)+1),0.0,0.0,0.0)
    f.write(line)

    for i in range(len(Ispin_info)):
        line = "{:>10.5f}\t{:>10.5f}\t{:>10.5f}\t{:>10}\n".format(\
        Ispin_info[i][0], Ispin_info[i][1], Ispin_info[i][2]            ,only_spin[i])
        f.write(line)
else :
    print("No atom\n")
    line = "{:>10.5f}\t{:>10.5f}\t{:>10.5f}\t{:>10.5f}\n".format(float(1),0.0,0.0,0.0)
    f.write(line)
f.close()

#4-2. Save the Defect info ####################################
print("\nNow Save the Defect information!-----------------------")
f= open(outputfile+'_Defect','w')

if (len(Defect_info_tmp) == 2):
    f.write("the vacancy site ({})\n".format(Vacancy))
    line= "{:>10.5f}\t{:>10.5f}\t{:>10.5f}\n".format(Defect_info_V[0], Defect_info_V[1], Defect_info_V[2])
    f.write(line)
    print("the vacancy site ({})".format(Vacancy))
    print("{:>10.5f}\t{:>10.5f}\t{:>10.5f}".format(Defect_info_V[0], Defect_info_V[1], Defect_info_V[2]))

f.write("the impurity site ({})\n".format(Impurity))
line= "{:>10.5f}\t{:>10.5f}\t{:>10.5f}\n".format(Defect_info_N[0], Defect_info_N[1], Defect_info_N[2])
print("the impurity site ({})".format(Impurity))
print("{:>10.5f}\t{:>10.5f}\t{:>10.5f}".format(Defect_info_N[0], Defect_info_N[1], Defect_info_N[2]))

f.write(line)

if (Opt_2D):
    f.write("\nCaution of using Bath size in 2D!\n")
    print("Caution of using Bath size in 2D!")
else:
    f.write("\nCaution of using Bath size!\n")
    print("Caution of using Bath size!")
 
f.write("sizeX : [{},{}] --> min Bath sizeX : {:>10.5f}\n".format(0,sizeX,min([sizeX-Defect_info_N[0],Defect_info_N[0]])))
f.write("sizeY : [{},{}] --> min Bath sizeY : {:>10.5f}\n".format(0,sizeY,min([sizeY-Defect_info_N[1],Defect_info_N[1]])))
f.write("sizeZ : [{},{}] --> min Bath sizeZ : {:>10.5f}\n".format(0,sizeZ,min([sizeZ-Defect_info_N[2],Defect_info_N[2]])))
print("sizeX : [{},{}] --> min Bath sizeX : {:>10.5f}".format(0,sizeX,min([sizeX-Defect_info_N[0],Defect_info_N[0]])))
print("sizeY : [{},{}] --> min Bath sizeY : {:>10.5f}".format(0,sizeY,min([sizeY-Defect_info_N[1],Defect_info_N[1]])))
print("sizeZ : [{},{}] --> min Bath sizeZ : {:>10.5f}".format(0,sizeZ,min([sizeZ-Defect_info_N[2],Defect_info_N[2]])))

 
if (Opt_2D):
    Bathlimit=[min([sizeX-Defect_info_N[0],Defect_info_N[0]]),min([sizeY-Defect_info_N[1],Defect_info_N[1]])]
else:
    Bathlimit=[min([sizeX-Defect_info_N[0],Defect_info_N[0]]),min([sizeY-Defect_info_N[1],Defect_info_N[1]]),min([sizeZ-Defect_info_N[2],Defect_info_N[2]])]
f.write("\nDo not jump the Bath size (limit : {:>10.5f})\n".format(min(Bathlimit)))
print("\nDo not jump the Bath size (limit : {:>10.5f})\n".format(min(Bathlimit)))

f.write("\nMade time : "+time.strftime('%c',time.localtime(time.time())))
f.close()

#4-3. Save the information info format : POSCAR (for test) ####################################
#test_on=True #to see the POSCAR file or structure
test_on=False #do not see the POSCAR

if (test_on):
    print("\nNow Save the POSCAR file !")
    f = open("Convert_2_POSCAR_form_"+outputfile+"_%d"%(q), 'w')
    print(only_spin2)
    
    ALL_namelist=[]
    print("ALL atomic_name")
    ALL_namelist.extend(atomic_name)
    ALL_namelist.extend(insert_atomic_name)
    print(ALL_namelist)
    Checking_num_Impu=ALL_namelist.count(Impurity)
    if(Vacancy != 'NONE'):
        Checking_num_Vacan=ALL_namelist.count(Vacancy)
    ALL_namelist=list(set(ALL_namelist))

    if Checking_num_Impu == 1:
        ALL_namelist.remove(Impurity)
    if(Vacancy != 'NONE'):
        if Checking_num_Vacan == 1:
            ALL_namelist.remove(Vacancy)

    print(ALL_namelist)
    
    f.write("ntyp = {}\n".format(len(ALL_namelist)))
    f.write("nat = {}\n".format((len(Ispin_info))))
    f.write("ATOMIC_SPECIES\n")
    for k in range(len(ALL_namelist)):
        f.write("{}\n".format(ALL_namelist[k]))
    
    f.write("CELL_PARAMETERS angstrom\n")
    f.write("{:>10.5f}\t0.00000\t0.00000\n".format(sizeX))
    f.write("0.00000\t{:>10.5f}\t0.00000\n".format(sizeY))
    f.write("0.00000\t0.00000\t{:>10.5f}\n".format(sizeZ))
    f.write("ATOMIC_POSITIONS (crystal)\n")

    if(len(only_spin2) != 0):

        only_spin2=only_spin2[:,3]
        Ispin_info=np.array(Ispin_info)
        Ispin_info=Ispin_info[:,:3].astype(np.float)

        p=re.compile("[^0-9]")
        for i in range(len(Ispin_info)):
            line = "{}\t{:>10.7f}\t{:>10.7f}\t{:>10.7f}\n".format(\
            "".join(p.findall(only_spin2[i])),Ispin_info[i][0]/sizeX, Ispin_info[i][1]/sizeY, Ispin_info[i][2]/sizeZ)

            f.write(line)
    else :
        print("No atom\n")
        line = "{:>10.5f}\t{:>10.5f}\t{:>10.5f}\t{:>10.5f}\n".format(float(1),0.0,0.0,0.0)
        f.write(line)
    f.close()

print("\nFinish to make the Bath!!------------------------------")
print("Now Making time : "+time.strftime('%c',time.localtime(time.time())))
