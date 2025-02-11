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
#2-1. ./BathGenerator_curve_insert.py [periodic inputfile] [Defect inputfile]  [insert Curvefile] [probability file] [size] [Impurity atom Name] [Vacancy atom Name] [outfile name] [inital number] [range number]
#
#--> only make the [size, size, size] cubic outfile using arbitray structure(POSCAR)
#--> if you don't want to make vaccancy atom, write the Vacancy atom Name = "NONE"
#
#
#2-2. ./BathGenerator_curve_insert.py [periodic inputfile] [Defect inputfile] [insert Curvefile] [probability file] [size] [Impurity atom Name] [Vacancy atom Name] [outfile name] [inital number] [range number] [2D option]
#--> [2D option] = '2D' or '2' (including "2")
#
#
###############################################################
#periodic inputfile     :   POSCAR which has bulk or pure structure
#Defect inputfile       :   POSCAR or csv file having the defect atom
#insert Curvefile       :   POSCAR which is curve of bulk or pure struture
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
        print("atom_data :",atom_data)
        print("atom_data_num :",atom_data_num)
        print("plz, check the atom species and number!!\n")
        sys.exit()
    elif len(unit_cell_parameter) != 3:
        print("unit_cell_parameter :",unit_cell_parameter)
        print("plz, check the cell parameter!!\n")
        sys.exit()
    elif len(temp_unit_cell_information) != sum(atom_data_num):
        print("temp_unit_cell_information :",temp_unit_cell_information)
        print("atom_data_num :",atom_data_num)
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


#2. solve the matrix to know the lattice index
############################
def MatSolv(latticeParameter,insert_latticeParameter,A,B,C):
#
#(cubic X, cubic Y, cubic Z) = (A, B, C)
# A = a1*(i-1) + b1*(j-1) + c1*(k-1) + a1' + b1' + c1'
# B = a2*(i-1) + b2*(j-1) + c2*(k-1) + a2' + b2' + c2'
# C = a3*(i-1) + b3*(j-1) + c3*(k-1) + a3' + b3' + c3'
#
#   latticeParameter(periodic file)
#   a1, a2, a3
#   b1, b2, b3
#   c1, c2, c3
#
#   insert_latticeParameter(inserted defect file)
#   a1', a2', a3'
#   b1', b2', b3'
#   c1', c2', c3'

    # a1, b1, c1       i-1    a1', b1', c1'     1     A
    # a2, b2, c2  dot  j-1  + a2', b2', c2' dot 1  =  B
    # a3, b3, c3       k-1    a3', b3', c3'     1     C
    #
    # [i j k] is # of a,b,c
    # [A B C] is length (unit : angstrom)
    #
    #--> knowing [i,j,k] ==> result

    Lattic=np.array(latticeParameter)
    LatticT=Lattic.T
    inv_Lat=lin.inv(LatticT)

    totA=A+(sum(latticeParameter[:][0]))-(sum(insert_latticeParameter[:][0]))
    totB=B+(sum(latticeParameter[:][1]))-(sum(insert_latticeParameter[:][1]))
    totC=C+(sum(latticeParameter[:][2]))-(sum(insert_latticeParameter[:][2]))

    Position=[[totA],[totB],[totC]]
    result=np.dot(inv_Lat, Position)

    return result

#3. Defect information -> remove the overlap one.
############################
def Defect_position(insert_unit_cell_information, Impurity, Vacancy,latticeParameter,insert_latticeParameter,sizeX,sizeY,sizeZ,Opt_2D,curve_lat,Defect_unitcell):

    #now the insert_unit_cell_information is str type due to atomic name
    #And the atomic name is 
    atomic_name=np.array(insert_unit_cell_information)
    atomic_name=atomic_name[:,0]


    #Find the mid lattice numbering 
    #Lattic=np.array(latticeParameter)
    #LatticT=Lattic.T
    #inv_Lat=lin.inv(LatticT)
    #Position=[[sizeX/2],[sizeY/2],[sizeZ/2]] #Temporary Defect position ,unit : angstrom
    #Mid_point=np.dot(inv_Lat, Position)
        
    MidA=Defect_unitcell[0]
    MidB=Defect_unitcell[1]
    MidC=Defect_unitcell[2]

    #in the insert curve code, we have to set the insert curve file.
    #it is difficult to decide the defect point unit cell
    #So, I just set mid unit cell point
    #if you have dissatisfaction on the mid unit cell point, 
    #Change the mid unit cell point method yourself. I don't have an idea another method

    #a_list = range(int(MidA-2)+1,math.ceil(MidA+2))
    #b_list = range(int(MidB-2)+1,math.ceil(MidB+2))
    #c_list = range(int(MidC-2)+1,math.ceil(MidC+2))
    a_list = [MidA]
    b_list = [MidB]
    c_list = [MidC]

    print("a_list(in Defect) : ",a_list)
    print("b_list(in Defect) : ",b_list)
    print("c_list(in Defect) : ",c_list)

    ref_length = lin.norm(np.array(insert_latticeParameter[2]))    #a3 length (angstrom)

    #Defect_p    = np.array([0.5,0.5,0.5]) #set ref position
    Defect_p    =np.array([sizeX/2,sizeY/2,sizeZ/2]) #ref position
    tmpDistance = ref_length+10 # angstrom
    n = np.array([0.0,0.0,0.0]) # to find close spin central
    Plus_lat=[0,0,ref_length] # add ref vacancy position (h_BN)

    #if the matrial is 2D, we need to check for the defect position in the 0-th layer
    p=re.compile("2")
    if (Opt_2D):
        #ref position
        Defect_p  =np.array([sizeX/2,sizeY/2,sizeZ/2]) 
        c_list=range(1)
        ref_lengthx=lin.norm(np.array(insert_latticeParameter[0]))
        ref_lengthy=lin.norm(np.array(insert_latticeParameter[1]))
        Plus_lat=[0,ref_length,0] # add ref vacancy position
        if (ref_lengthx < ref_lengthy):
            tmpDistance = ref_lengthy # angstrom
        else:
            tmpDistance = ref_lengthx # angstrom

#    print("NV position candidate list :\n")
    for a in a_list:
        for b in b_list:
            for c in c_list:
                for p0 in range(0,len(insert_unit_cell_information)):
                    if atomic_name[p0] == Impurity:
                        Data_tmp=[a,b,c,p0]
                        tmpN=RealPosition(latticeParameter,insert_latticeParameter,insert_unit_cell_information,Data_tmp,Defect_unitcell,curve_lat)

                        distance = dist(tmpN,Defect_p)
                        if distance <= tmpDistance:
                            print(1)
                            tmpDistance = distance 
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
        checkDistance = ref_length #angstrom
    
        print('origin vacancy       : ', v_tmp)

    #the vacancy is in the same lattice
        a_list=[N_unit_data[0]]
        b_list=[N_unit_data[1]]
        c_list=[N_unit_data[2]]
    # test the v site is real
        for a in a_list:
            for b in b_list:
                for c in c_list:
                    for p0 in range(0,len(insert_unit_cell_information)):
                        if (atomic_name[p0] == Vacancy and p0 != N_unit_data[3]) :
                            Data_tmp1=[a,b,c,p0]
                            tmpV=RealPosition(latticeParameter,insert_latticeParameter,insert_unit_cell_information,Data_tmp1,Defect_unitcell,curve_lat)

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

#4. read the nucelar spin probability of isotopes nuclear
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
    
    f.close()
        
    return pro_name, pro_spin, pro_num


#5.real position of atom 
############################
def RealPosition(latticeParameter,insert_latticeParameter,unit_cell_information,data_set,Defect_unitcell,curve_lat):
    LengthX     = np.array(latticeParameter[0])
    LengthY     = np.array(latticeParameter[1])
    LengthZ     = np.array(latticeParameter[2])

    LengthXD     = np.array(insert_latticeParameter[0])
    LengthYD     = np.array(insert_latticeParameter[1])
    LengthZD     = np.array(insert_latticeParameter[2])

    A=data_set[0]
    B=data_set[1]
    C=data_set[2]
    p0=data_set[3]


    #(real site)
    #curve_lat == 3, 4, 5
    #Real = (unit_cell_information[p0][1]+A)*LengthX + (unit_cell_information[p0][2]+B)*LengthY + (unit_cell_information[p0][3]+C)*LengthZ 
    RealX = (unit_cell_information[p0][1]+A)*LengthX 
    RealY = (unit_cell_information[p0][2]+B)*LengthY 
    RealZ = (unit_cell_information[p0][3]+C)*LengthZ 

    if curve_lat == 0 or curve_lat == 1 or curve_lat ==2:
        if (data_set[curve_lat] > Defect_unitcell[curve_lat]):
            if(curve_lat == 0):
                #Real = (unit_cell_information[p0][1]+A-1)*LengthX + (unit_cell_information[p0][2]+B)*LengthY + (unit_cell_information[p0][3]+C)*LengthZ + LengthXD
                RealX = (unit_cell_information[p0][1]+A-1)*LengthX + LengthXD
            elif(curve_lat == 1):
                #Real = (unit_cell_information[p0][1]+A)*LengthX + (unit_cell_information[p0][2]+B-1)*LengthY + (unit_cell_information[p0][3]+C)*LengthZ + LengthYD
                RealY = (unit_cell_information[p0][2]+B-1)*LengthY + LengthYD
            elif(curve_lat == 2):
                #Real = (unit_cell_information[p0][1]+A)*LengthX + (unit_cell_information[p0][2]+B)*LengthY + (unit_cell_information[p0][3]+C-1)*LengthZ + LengthZD
                RealZ = (unit_cell_information[p0][3]+C-1)*LengthZ + LengthZD

        elif (data_set[curve_lat] == Defect_unitcell[curve_lat]):
            if(curve_lat == 0):
                #Real = (A)*LengthX + (unit_cell_information[p0][2]+B)*LengthY + (unit_cell_information[p0][3]+C)*LengthZ + LengthXD*unit_cell_information[p0][1]
                RealX = (A)*LengthX + LengthXD*unit_cell_information[p0][1] 
            elif(curve_lat == 1):
                #Real = (unit_cell_information[p0][1]+A)*LengthX + (B)*LengthY + (unit_cell_information[p0][3]+C)*LengthZ + LengthYD*unit_cell_information[p0][2]
                RealY = (B)*LengthY + LengthYD*unit_cell_information[p0][2] 
            elif(curve_lat == 2):
                #Real = (unit_cell_information[p0][1]+A)*LengthX + (unit_cell_information[p0][2]+B)*LengthY + (C)*LengthZ + LengthZD*unit_cell_information[p0][3]
                RealZ = (C)*LengthZ + LengthZD*unit_cell_information[p0][3] 

    Real = RealX + RealY + RealZ

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

if (len(sys.argv) != 11 and len(sys.argv) != 12 ):
    print("your option is poor!!")
    print("[1] : Periodic File, [2] : insert defect file name, [3] : insert Curve file name,  [4] : probability File, [5] : bath size(cubic), [6] : Impurity atom name, [7] : Vacancy atom name, [8] : save file name, [9] : inital file number, [10] : last file number, [11] : 2D option ")

Periodicfile    = sys.argv[1]
insertDefectfile= sys.argv[2]
insertCurvefile = sys.argv[3]

concTempfile    = sys.argv[4]
size            = int(sys.argv[5])  #angstrom
                                    #make the lattice as cubic.
Impurity        = sys.argv[6]
Vacancy         = sys.argv[7]
outputfile      = sys.argv[8]
numOfConfigIni  = int(sys.argv[9])  #inital numbering
numOfConfig     = int(sys.argv[10])  #range output data file

if (len(sys.argv) == 12):
    print(len(sys.argv))
    Option2D       = sys.argv[11]
    Opt_2D=ReadOption(Option2D)
    print("2D_axis : ",Opt_2D)
else :
    Opt_2D=False
print("2D_axis option : ",Opt_2D)
    

#2. Other parameter ##################################
print('\n----------- Inputfile and Setting number!! ----------\n')
con_unit            = 0.01 #due to [unit : %] (1% = 1000000ppm)
                           #convert [% unit] to numbering
#about Periodic inputfile
print('\nRead bulk file!')
latticeParameter, atom_data, atom_data_num, unit_cell_information = unitCellInfo(Periodicfile)

#Insert Defect file 
print('\nRead insert file!')
insert_latticeParameter, insert_atom_data, insert_atom_data_num, insert_unit_cell_information = unitCellInfo(insertDefectfile)
#print("insert_unit_cell_information : ", insert_unit_cell_information)

#Insert Curve file 
print('\nRead insert Curve file!')
insertCurve_latticeParameter, insertCurve_atom_data, insertCurve_atom_data_num, insertCurve_unit_cell_information = unitCellInfo(insertCurvefile)
#print("insertCurve_unit_cell_information : ", insertCurve_unit_cell_information)

print("latticeParameter")
print(latticeParameter)
print("insert_latticeParameter")
print(insert_latticeParameter)
print("insertCurve_latticeParameter")
print(insertCurve_latticeParameter)

#Check the lattice parameter is correct or what is wrong!
#Must, the Defect file and Curve file have to be same in lattice parameter
if (insert_latticeParameter != insertCurve_latticeParameter):
    print("Maybe, The lattice parameter of defect file and insert Curve file wrong!")
    print("Check the defect file and Curve file to be same in lattice parameter!")
    sys.exit()

#[0] : a-axis & [1] : b-axis & [2] : c-axis
#[0] : x, [1] : y, [2] : z, [3] : xy, [4] : xz, [5] : yz
check_lat_list=[]
for i in range(2):
    if (latticeParameter[i] != insert_latticeParameter[i]):
        print("the difference is :", i)
        check_lat_list.append(i)
if len(check_lat_list) == 0:
    print("you can choose the curve_lat 0~5")
    print("[0] : x, [1] : y, [2] : z, [3] : xy, [4] : xz, [5] : yz")
    print("--> [0] : x-plane, [1] : y-plane, [2] : z-plane,")
    print("--> [3] : along z-axis, [4] : along y-axis, [5] : along x-axis")
    curve_lat=3  #3 : xy --> curve file stack along z-axis
    print("Now you choose insert curve file format :",curve_lat)
elif len(check_lat_list) == 1:
    curve_lat = check_lat_list[0]
    print("Now you are only able to choose insert curve file format :",curve_lat)
    print("if you want to choose curve_lat=3~5, uniform the lattice parameter!!")
else:
    print("The defect file and insert file must have one different lattice parameter, or the same")
    print("Because we have to satisfy the periodic condition")

    for i in range(len(check_lat_list)):
        print("latticeParameter{0}".format([check_lat_list[i]]))
        print(latticeParameter[check_lat_list[i]])
        print("insert_latticeParameter{0}".format([check_lat_list[i]]))
        print(insert_latticeParameter[check_lat_list[i]])
        print()
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
            count_lat=MatSolv(latticeParameter,insert_latticeParameter,i,j,k)
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

#set the defect unit cell using MatSolv
tmpcal=(MatSolv(latticeParameter,insert_latticeParameter,sizeX,sizeY,sizeZ))/2
Defect_unitcell=[int(tmpcal[0]),int(tmpcal[1]),int(tmpcal[2])]
print("Defect_unitcell")
print(Defect_unitcell)


#print('\nunitcell info : \n',unit_cell_information)
print('\nlatticeParameter(periodic lat)       : \n',latticeParameter)

#Find the Defect center spin position
print('\n----------- Defect site finding!! ----------\n')
#Defect_info_tmp, Defect_data = Defect_position(unit_cell_information, Impurity, Vacancy,latticeParameter,sizeX,sizeY,sizeZ,Opt_2D)
print("We find the Defect position in the",insertDefectfile)
Defect_info_tmp, Defect_data = Defect_position(insert_unit_cell_information, Impurity, Vacancy,latticeParameter,insert_latticeParameter,sizeX,sizeY,sizeZ,Opt_2D,curve_lat,Defect_unitcell)
print('Length of Defect_information : ', len(Defect_info_tmp))

N_info_tmp = Defect_info_tmp[0]
N_data_tmp = Defect_data[0]
print('\nN real site        : ',N_info_tmp)
print('\nN data             : ',N_data_tmp)

#Impurity site site (real site)
Defect_info_N=RealPosition(latticeParameter,insert_latticeParameter,insert_unit_cell_information,N_data_tmp,Defect_unitcell,curve_lat) 

print('Defect of impurity : ',Defect_info_N)
Defect_site_tmp=[]
Defect_site_tmp.append(N_info_tmp.tolist())


if (len(Defect_info_tmp) == 2):
    V_info_tmp = Defect_info_tmp[1]
    V_data_tmp = Defect_data[1]
    print('\nV real site        : ',V_info_tmp)
    print('V data             : ',V_data_tmp)

    #Vacancy site (real site)
    Defect_info_V=RealPosition(latticeParameter,insert_latticeParameter,insert_unit_cell_information,V_data_tmp,Defect_unitcell,curve_lat)  

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
    elif(iso_name[0] not in atom_data) and (iso_name[0] not in insert_atom_data) and (iso_name[0] not in insertCurve_atom_data):
        if(iso_name[0] is not Impurity) or (iso_name[0] is not Vacancy):
            print("\n####################################################\n")
            print("The atom name of isotope is not in the unit cell\n")
            print("plz, check the isotope file(probablity file)\n")
            print("other atoms?\n")
            print("Checking iso name :", iso_name[0])
            print("\n####################################################\n")
            sys.exit()
        else:
            print("Is this the Vacancy or Impurity?") 
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
        elif iso_name[0] in insertCurve_atom_data:
            index_num=insertCurve_atom_data.index(iso_name[0])
            isotope_data_num.append(insertCurve_atom_data_num[index_num])
            isotope_data_name.append(insertCurve_atom_data[index_num])

print('Isotope data Name : ', isotope_data_name)
print('Isotope data Number : ', isotope_data_num)


#3. Random generator spin within supercell size
print('\n----------- Random generator!! ----------\n')

#atomic name
atomic_name=np.array(unit_cell_information)
atomic_name=atomic_name[:,0]

insert_atomic_name=np.array(insert_unit_cell_information)
insert_atomic_name=insert_atomic_name[:,0]

insertCurve_atomic_name=np.array(insertCurve_unit_cell_information)
insertCurve_atomic_name=insertCurve_atomic_name[:,0]

loopX=range(unitCellnumber[0][0], unitCellnumber[0][1]+1)
loopY=range(unitCellnumber[1][0], unitCellnumber[1][1]+1)
loopZ=range(unitCellnumber[2][0], unitCellnumber[2][1]+1)

if (Opt_2D):
    loopZ=range(1)
#checking for configure atom
Checking_atom=0

if curve_lat == 0:
    curve_lat1=0
    curve_lat2=0
elif curve_lat == 1:
    curve_lat1=1
    curve_lat2=1
elif curve_lat == 2:
    curve_lat1=2
    curve_lat2=2
elif curve_lat == 3:
    curve_lat1=0
    curve_lat2=1
elif curve_lat == 4:
    curve_lat1=0
    curve_lat2=2
elif curve_lat == 5:
    curve_lat1=1
    curve_lat2=2

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
                if curve_lat == 0: #x
                    curve_unit1 = i
                    curve_unit2 = i

                elif curve_lat == 1: #y
                    curve_unit1 = j
                    curve_unit2 = j

                elif curve_lat == 2: #z
                    curve_unit1 =k
                    curve_unit2 =k

                elif curve_lat == 3: #xy
                    curve_unit1 =i
                    curve_unit2 =j

                elif curve_lat == 4: #xz
                    curve_unit1 =i
                    curve_unit2 =k

                elif curve_lat == 5: #yz
                    curve_unit1 =j
                    curve_unit2 =k

                if Onprint:print("SuperCell[%d,%d,%d], Ispin_info_temp \n:"%(i,j,k))
                #About the defect atom unit cell
                if (i == N_data_tmp[0]) and (j==N_data_tmp[1]) and (k==N_data_tmp[2]): #the defect file atom

                    for n in range(len(insert_unit_cell_information)):
                        Ispin_info_temp = []
                        Ispin_data_temp=[i,j,k,n]
                        Ispin_info_temp=RealPosition(latticeParameter,insert_latticeParameter,insert_unit_cell_information,Ispin_data_temp,Defect_unitcell,curve_lat)
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

                elif (curve_unit1 == Defect_unitcell[curve_lat1]) and (curve_unit2 == Defect_unitcell[curve_lat2]) : #the curve file atom
                    for n in range(len(insertCurve_unit_cell_information)):
                        Ispin_info_temp = []
                        Ispin_data_temp=[i,j,k,n]
                        Ispin_info_temp=RealPosition(latticeParameter,insert_latticeParameter,insertCurve_unit_cell_information,Ispin_data_temp,Defect_unitcell,curve_lat)
                        if Onprint:{print("Ispin data info : ",Ispin_data_temp),
                                    print("atomic name     : ",insertCurve_atomic_name[n]),
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
                            p=re.compile(insertCurve_atomic_name[n])
    
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


                else:
                    for n in range(len(unit_cell_information)):
                        Ispin_info_temp = []
                        Ispin_data_temp=[i,j,k,n]
                        Ispin_info_temp=RealPosition(latticeParameter,insert_latticeParameter,unit_cell_information,Ispin_data_temp,Defect_unitcell,curve_lat)
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
                                                print("pro_name[n] :",pro_name[Iso_count[SetIso[0]]]) }
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
    Ispin_info=Ispin_info[:,:3].astype(np.float64)

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
    ALL_namelist.extend(insertCurve_atomic_name)
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
