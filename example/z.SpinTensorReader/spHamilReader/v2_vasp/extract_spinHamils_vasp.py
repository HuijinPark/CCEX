import numpy as np
import sys
from spHamilReader.base.rwfiles import readfile
from spHamilReader.v1_vasp.read_vasp import *
from spHamilReader.base.printatm import *
from spHamilReader.base.ioPoscar import read_POSCAR
from spHamilReader.base.dictproc import print_dict
from spHamilReader.spin_database.loadisotopes import isotopeData

def readtensor_3x3(linedata):
    tensor = list(map(float,linedata))
    return    [tensor[0],tensor[3],tensor[4],\
             tensor[3],tensor[1],tensor[5],\
             tensor[4],tensor[5],tensor[2]]

def readtensor_diag(linedata):
    tensor = list(map(float,linedata))
    return    [tensor[0],tensor[1],tensor[2]]

def readtensor_elements(linedata,elem):
    tensor = list(map(float,linedata))
    return    [tensor[i] for i in range(elem)]

def readtensor_Aiso(linedata,addhfcore=True):
    """
    frozencore : A_1c in VASP
    fermiconta : A_tot in VASP
        See https://www.vasp.at/wiki/index.php/LHYPERFINE
        You should add the core electronic contribution
    """
    coreeffect=0
    if addhfcore : coreeffect = float(linedata[3])
    fermiconta = float(linedata[4])
    return fermiconta + coreeffect

def printform(elements,contents,_list):
    strform="\t{:>7}"+"{:>10}"*elements + "\n\n"
    valform="\t{:>7}"+"{:>10}"*elements + "\n" 
    message = ""
    message += strform.format("idx",*contents)

    for i in range(len(_list)):
        for j in range(len(_list[i])):
            _list[i] = list(_list[i])
            try:
                _list[i][j] = "{:>10.3f}".format(_list[i][j])
            except:
                _list[i][j] = "{:>10}".format(_list[i][j])

    if len(_list) > 10:
        for j in range(3):
            message += valform.format(f"({j+1})",*_list[j])
        message += "\t{:>7}\n".format(":")
        message += valform.format(f"({len(_list)-2+1})",*_list[len(_list)-2])
        message += valform.format(f"({len(_list)-1+1})",*_list[len(_list)-1])
    else:
        for j in range(len(_list)):
            message += valform.format(f"({j+1})",*_list[j])
    return message

def extQdata(atomicdata,dat,print_io=False,startingline=0):

    keys={'EFG'         :[],    # ion# x [Vxx, Vxy, Vxz, Vyx .., Vzz] (V/A^2)
          'EFG_diag'    :[],    # ion# x [Vxx(D), Vyy(D), Vzz(D)] (V/A^2) 
          'Eigenvectors':[],    # ion# x [vxx, vxy, vxz, vyx .., vzz]
          'Cq'          :[],    # ion# x Cq (MHz)
          'eta'         :[],    # ion# x eta 
          'q'           :[],    # ion# x Q (mb)
          'Cq,eta,q'    :[]}
    
    ###############################################
    read=False
    line=startingline
    while(line<len(dat)):
        ###############################################
        # find EFG tensor
        if checkLine_all(dat[line]," Electric field gradients (V/A^2)"):
            readtensor=readtensor_3x3
            key='EFG'
            read=True
        # find EFG tensor diag
        if checkLine_all(dat[line]," Electric field gradients after diagonalization (V/A^2)"):
            readtensor=readtensor_diag
            key='EFG_diag'
            read=True
        ## find Eigenvector (How to..?)
        #if checkLine_all(dat[line]," Eigenvectors"):
        #    readtensor=readtensor_eigen
        #    key='Eigenvectors'
        #    read=True
        # find Cq, eta,q
        if checkLine_all(dat[line],"NMR quadrupolar parameters"):
            readtensor=readtensor_diag
            key='Cq,eta,q'
            read=True
        ###############################################
        if read == True:
            nion=0
            while(nion<atomicdata["total_atom_num"]):
                try:
                    keys[key].append(readtensor(dat[line][1:]))
                    nion +=1
                except: pass;
                line+=1
            read=False
        line+=1
        ###############################################
        #break while when you got all materials
        if keys['Cq,eta,q'] != []:
            break;
    ###############################################

    tmp = np.array(keys['Cq,eta,q'])
    keys['Cq']  = tmp[:,0]
    keys['eta'] = tmp[:,1]
    keys['q']   = tmp[:,2]
    del keys['Cq,eta,q']

    if print_io == True:
        for k,v in keys.items():    
            print("-"*40)
            if k=='EFG' or k == 'EFG_diag':
                print(k,'(V/A^2')
            if k=='Cq':
                print(k,'= e*Q*V_zz/h (MHz)')
            if k=='eta':
                print(k,'=(V_yy-V_xx)/V_zz')
            if k=='q':
                print(k,'(mb = 10^-31)')
            print("-"*40)
            for j in range(len(v)):
                print(j, v[j])
            print("-"*40)

    return (line,keys)

def extAdata(atomicdata,dat,print_io=False,startingline=0,Ngyro=None,_addhfcore=True):

    keys={'Aiso'    :[],    # ion# x Fermi contact(MHz) 
          'Adip'    :[],    # ion# x [Axx, Axy, Axz, Ayx .., Azz] (MHz)
          'Adiag'   :[],    # ion# x [Axx(D), Ayy(D), Azz(D)] (MHz) for total A
          'givenr'  :{}}    # ion# x given gamma (MHz/T) in incar 

    ###############################################
    read=False
    line=startingline
    while(line<len(dat)):
        ###############################################
        # find Fermi contact 
        if checkLine_all(dat[line],\
            " Fermi contact (isotropic) hyperfine coupling parameter (MHz)"):
            readtensor=readtensor_Aiso
            key='Aiso'
            read=True
        # find A dipolar term 
        if checkLine_all(dat[line],\
            "Dipolar hyperfine coupling parameters (MHz)"):
            readtensor=readtensor_3x3
            key='Adip'
            read=True
        # find diagonalized term
        if checkLine_all(dat[line],\
            "Total hyperfine coupling parameters after diagonalization (MHz)"):
            readtensor=readtensor_diag
            key='Adiag'
            read=True
        ###############################################
        if read == True:
            nion=0
            while(nion<atomicdata["total_atom_num"]):
                try:
                    if key=="Aiso":
                        keys[key].append(readtensor(dat[line][1:],addhfcore=_addhfcore))
                    else:
                        keys[key].append(readtensor(dat[line][1:]))
                    nion +=1
                except: pass;
                line+=1
            read=False
        line+=1
        ###############################################
        #break while when you got all materials
        if keys['Adiag'] != []:
            break;
    ###############################################

    if Ngyro!=None:
        Ngyro = list(map(float,Ngyro.strip().split()))
        for i, atm in enumerate(atomicdata['atom_data']):
            keys['givenr'][atm] = Ngyro[i]
    else:
        for i, atm in enumerate(atomicdata['atom_data']):
            keys['givenr'][atm] = 1.0

    ###############################################
    if print_io == True:
        for k,v in keys.items():    
            print("-"*40)
            print(k,"(MHz)")
            print("-"*40)
            if k!='givenr':
                for j in range(len(v)):
                    print(j, v[j])
            else:
                print_dict(v)
            print("-"*40)

    return (line,keys)

def extDdata(dat,print_io=False):

    keys={'ZFS'        :[],    # np.array([Dxx, Dxy, Dxz, Dyx .., Dzz]) (MHz)
          'D'        :[],    # 3/2*Dzz (MHz)
          'E'        :[],    # (Dxx - Dyy)/2 (MHz)
          'D_diag'    :[],    # [Dxx(D),Dyy(D),Dzz(D)] (MHz) 
          'D_eigvec':[]}    # [[vec_x, vec_y, vec_z],[],[]]
    
    ###############################################
    dat=dat[-150:]    
    line=0
    while(line<len(dat)):
        ###############################################
        if checkLine_all(dat[line],\
             "Spin-spin contribution to zero-field splitting tensor (MHz)"):
            nion=0
            while(nion<1):
                try:
                    keys['ZFS']=readtensor_3x3(dat[line])
                    nion +=1
                    break;
                except: pass;
                line+=1
            keys['ZFS']=np.array(keys['ZFS']).reshape(3,3)

        if checkLine_all(dat[line],\
             "after diagonalization"):
            nion=0
            while(nion<3):
                try:
                    tmp = readtensor_elements(dat[line],4)
                    keys['D_diag'].append(tmp[0])
                    keys['D_eigvec'].append(tmp[1:])
                    nion +=1
                except: pass;
                line+=1
            keys['D'] = (3.0/2.0)*keys['D_diag'][2]
            keys['E'] = (keys['D_diag'][0] - keys['D_diag'][1])/2.0
            
        ###############################################
        line+=1
        ###############################################
        #break while when you got all materials
        if keys['D_diag'] != []:
            break;
    ###############################################

    if print_io ==True:
        for k,v in keys.items():    
            print("-"*40)
            if k=='ZFS' or k=='D_diag':
                print(k,"(MHz)")
            if k=='D':
                print(k,"=3/2*D_diag_zz(MHz)")
            if k=='E':
                print(k,"=(D_diag_xx-D_diag_yy)/2(MHz)")
            if k=='D_eigvec':
                print(k)
            print("-"*40)
            print(v)
            print("-"*40)
    return keys


class extSpinHamilTensors:

    """
        Args:
            foutcar(str)    : outcar file
                              this is the main file that include spin Hamiltonian tensor
                              We will use last iteraction data
            fincar(str,opt)    : incar file (optional)
                              this file is used to check if you run DFT calculation as spin Hamiltonian option is True
                              
            readA(bool,opt)    : True or False to read A tensor or not 
            readQ(bool,opt)    : True or False to read Q tensor or not
            readD(bool,opt)    : True or False to read D tensor or not
                                above three option decide which tensor you will obtain
                                But if you turn the "fincar" on, fincar would be priority
        Return:
            * incar(dict)  : INCAR data from OUTCAR 
            * outcar(list) : OUTCAR data from OUTCAR 
            * atomicdata(dict) : Simple poscar from OUTCAR
            * Qdata(dict)   
            * Adata(dict)      
            * Ddata(dict)
            * isodata(dict)

            A,Q,D tensors(Dict)
            
            Keys:

            # CAUTION : ALL Materials are Raw data of outcar

            For Q tensor : Qdata 
                'EFG'           :[]        # ion# x [Vxx, Vxy, Vxz, Vyx .., Vzz] (V/A^2)
                'EFG_diag'      :[]        # ion# x [Vxx(D), Vyy(D), Vzz(D)] (V/A^2) 
                'Cq'            :np.arr    # ion# x Cq (MHz)
                'eta'           :np.arr    # ion# x eta 
                'q'             :np.arr    # ion# x Q (mb)
                ## None ! 'Eigenvectors' :[] # ion# x [vxx, vxy, vxz, vyx .., vzz] 

            For A tensor : Adata
                'Aiso'          :[]        # ion# x Fermi contact(MHz) 
                'Adip'          :[]        # ion# x [Axx, Axy, Axz, Ayx .., Azz] (MHz)
                'Adiag'         :[]        # ion# x [Axx(D), Ayy(D), Azz(D)] (MHz) for total A
                'givenr'        :[]        # ion# x given gamma (MHz/T) in incar 

            For D tensor : Ddata
                'ZFS'           :np.arr    # np.array([Dxx, Dxy, Dxz, Dyx .., Dzz]) (MHz)
                  'D'           :float     # 3/2*Dzz (MHz)
                'E'             :float     # (Dxx - Dyy)/2 (MHz)
                'D_diag'        :[]        # [Dxx(D),Dyy(D),Dzz(D)] (MHz) 
                'D_eigvec'      :[]        # [[vec_x, vec_y, vec_z],[],[]]

            Isotope information : isodata
                key(str) : atom name (ex. B)
                value(myspindict) : isotope information for each atom (ex. Isodata['B'] =  {'10B': (.. , .., ) })
                values's key : isotope name (ex. '10B')
                values's value : isotope information  (ex. (spin, gyro, q, concent ) ) 
                ==> summary : Isodata['B'] =  {'10B': (3.0, 2.8747, 0.08450, 0.1990) }
                isodata['B']['10B'].s          : spin         (no unit)
                isodata['B']['11B'].gyro      : gyro..     (rad/ms/G)
                isodata['N']['14N'].q          : eQ         (barn = 10^-28 m^2)
                isodata['N']['15N'].conc      : concent.. (%/100)
            
    """


    def __init__(self,foutcar,print_io=False,addhfcore=True):
        self.outcar, self.atomicdata, self.incar = readOUTCAR(foutcar,readincar=True)

        self.Qdata={}
        self.Adata={}
        self.Ddata={}
        self.isodata={}

        print("\tRead INCAR part in OUTCAR ...")

        try:
            self.readA = self.incar['LHYPERFINE']
        except: 
            self.readA = '.FALSE.'
        try:
            self.readQ = self.incar['LEFG']
        except:
            self.readQ = '.FALSE.'
        try:
            self.readD = self.incar['LDMATRIX']
        except:
            self.readD = '.FALSE.'

        print(f"\tRead A tensor : {self.readA}")
        print(f"\tRead Q tensor : {self.readQ}")
        print(f"\tRead D tensor : {self.readD}")

        try:
            self.NGYROMAG = self.incar['NGYROMAG']
            print(f"\tNGYROMAG : {self.NGYROMAG}")
            #sys.exit("\tERROR, NGYRO term need to be considered in a code")
        except:
            self.NGYROMAG = None
            print(f"\tNGYROMAG : NONE ")

        try:
            self.QUAD_EFG = self.incar['QUAD_EFG']
            print(f"\tQUAD_EFG : {self.QUAD_EFG}")
            sys.exit("\tERROR, QUAD_EFG term need to be considered in a code")
        except:
            self.QUAD_EFG = None
            print(f"\tQUAD_EFG : NONE ")


        print("\tRead POSCAR(simple) part in OUTCAR ...")
        print(f"\t atoms # = {self.atomicdata['total_atom_num']}")
        print(f"\t atomic species = {self.atomicdata['atom_data']}")

        line=0
        # OUTCAR always in order Q-A-D
        if 'T' in self.readQ :
            print("\tRead Quadrupole tensor part in OUTCAR ...")
            line,self.Qdata= extQdata(self.atomicdata,self.outcar,print_io=print_io)
        if 'T' in self.readA :
            print("\tRead Hyperfine tensor part in OUTCAR ...")
            line,self.Adata = extAdata(self.atomicdata,self.outcar,print_io=print_io,_addhfcore=addhfcore,\
            startingline=line,Ngyro=self.NGYROMAG)
        if 'T' in self.readD :
            print("\tRead ZFS tensor part in OUTCAR ...")
            self.Ddata = extDdata(self.outcar,print_io=print_io)

        c=0
        print("\t",self.atomicdata['atom_data'])
        #print(isotopeData[self.atomicdata['atom_data'][0]])
        for v in self.atomicdata['atom_data']:
            print("\tatm:",v)
            print("\tisotopeData",isotopeData(v))
            try:
                self.isodata[v] = isotopeData(v) 
                c+=1
            except:
                for i in range(c):
                    print("\t",i,self.atomicdata['atom_data'][i])
                    if self.atomicdata['atom_data'][i] in v:
                        self.isodata[v] = self.isodata[self.atomicdata['atom_data'][i]]
                        break;
            print("\tiso",v,self.isodata[v])
        
    def __repr__(self):
        message = "\n"
        message += "\t"+"-"*40 + "\n"
        message += f"\t{type(self).__name__} : \n"
        message += "\t"+"-"*40 + "\n"

        message += self.prt_atomicdata()

        if 'T' in self.readQ:
            message += self.prt_readQ()
        if 'T' in self.readA:
            message += self.prt_readA()
        if 'T' in self.readD:
            message += self.prt_readD()
        
        message += self.prt_isotopedata()

        return message

    def prt_readQ(self):
        ##################################
        EFG=self.Qdata['EFG']
        EFG_diag=self.Qdata['EFG_diag']
        Cq=self.Qdata['Cq']
        eta=self.Qdata['eta']
        q=self.Qdata['q']
        ###################################
        message = "" 
        message += "\n\tQuadrupole tensor (KEY : Qdata)\n"
        message += "\t"+"-"*30 + "\n"
        message += "\n\t*EFG : EFG tensor (V/A^2) \n\n"
        message += printform(9,[ "Vxx","Vxy","Vxz",\
                                 "Vyx","Vyy","Vyz",\
                                 "Vzx","Vzy","Vzz"],EFG)
        ###################################
        message += "\n\t*EFG_diag : Diagonalized EFG tensor (V/A^2) \n\n"
        message += printform(3,["Vxx(D)","Vyy(D)","Vzz(D)"],EFG_diag)
        ###################################
        message += "\n\tNMR quadrupolar parameters\n\n"
        message += "\t*Cq  : quadrupolar parameter    Cq=e*Q*V_zz/h (MHz) \n"
        message += "\t*eta : asymmetry parameters     (V_yy - V_xx)/ V_zz \n"
        message += "\t*Q   : nuclear electric quadrupole moment in mb (millibarn = 10^-31 m^2) \n\n"
        message += printform(3,["Cq","eta","q"],list(zip(Cq,eta,q)))
        message += "\t"+"-"*30 + "\n"
        ###################################
        return message

    def prt_readA(self):
        Aiso = self.Adata['Aiso']
        Adip = self.Adata['Adip']
        Adiag = self.Adata['Adiag']
        ###################################
        message = "" 
        message += "\n\tHyperfine tensor (KEY : Adata)\n"
        message += "\t"+"-"*30 + "\n"

        message += "\n\tHyperfine tensor \n"
        message += "\n\t*Aiso : Fermi contact term (MHz) \n"
        message += "\t*Adip : Dipolar coupling term (MHz) \n\n"
        message += printform(10,[ "Aiso","Axx","Axy","Axz",\
                                           "Ayx","Ayy","Ayz",\
                                           "Azx","Azy","Azz"],\
                                  [[fc]+dip for fc,dip in zip(Aiso,Adip)])
        ###################################
        message += "\n\t*Adiag : Diagonalized total Hyperfine tensor (MHz) \n\n"
        message += printform(3,[ "Axx(D)","Ayy(D)","Azz(D)"],Adiag)
        message += "\t"+"-"*30 + "\n"
        ###################################
        return message

    def prt_readD(self):
        message = "" 
        message += "\n\tZero Field splitting tensor (KEY : Ddata)\n"
        message += "\t"+"-"*30 + "\n"

        message += "\n\t*ZFS : ZFS tensor (MHz) \n\n"
        message += "\t{:>10.3f}{:>10.3f}{:>10.3f}\n".format(self.Ddata['ZFS'][0,0],self.Ddata['ZFS'][0,1],self.Ddata['ZFS'][0,2])
        message += "\t{:>10.3f}{:>10.3f}{:>10.3f}\n".format(self.Ddata['ZFS'][1,0],self.Ddata['ZFS'][1,1],self.Ddata['ZFS'][1,2])
        message += "\t{:>10.3f}{:>10.3f}{:>10.3f}\n".format(self.Ddata['ZFS'][2,0],self.Ddata['ZFS'][2,1],self.Ddata['ZFS'][2,2])
            
        message += "\n\tDiagonalized ZFS tensor \n"
        message += "\t*D_diag : Diagonalized ZFS tensor (MHz) \n"
        message += "\t*D_eigvec : Eigenvectors for diagonalized ZFS tensor \n\n"
        message += printform(4,[ "Dii(D)", 'eigvec_x', 'eigvec_y', 'eigvec_z'],\
            [[d]+v for d,v in zip(self.Ddata['D_diag'],self.Ddata['D_eigvec'])])

        message += "\n\tZFS parameters (MHz) \n"
        message += "\t*D : axial component D = 3/2*D_zz (MHz) \n"
        message += "\t*E : transverse component E = (D_xx-D_yy)/2 (MHz) \n\n"
        message += "\t{:>7}{:>10.3f}".format('D = ',self.Ddata['D']) + "\n"
        message += "\t{:>7}{:>10.3f}".format('E = ',self.Ddata['E']) + "\n"
        message += "\t"+"-"*30 + "\n"
        return message

    def prt_atomicdata(self):
        message = "" 
        message += "\n\tAtomic data from OUTCAR (KEY : atomicdata)\n"
        message += "\t"+"-"*30 + "\n"

        message += f"\t *atom_data : {self.atomicdata['atom_data']}\n" 
        message += f"\t *total_atom_num : {self.atomicdata['total_atom_num']}\n"
        message += f"\t *unit_cell_information :  \n\n"
        message += printform(3,["frac_x","frac_y","frac_z"],\
                                self.atomicdata["unit_cell_information"])
        return message

    def prt_isotopedata(self):
        message = "" 
        message += "\n\tIsotopic data (KEY : isodata)\n"
        message += "\t"+"-"*30 + "\n\n"

        message += "\t Atom : { Isotopes : ( S, gyro, eQ, conc) }\n\n"

        message += "\texample) \n"

        message += "\t isodata['B']['10B'].s          : spin          (no unit)\n"
        message += "\t isodata['B']['11B'].gyro       : gyro          (rad/ms/G)\n"
        message += "\t isodata['N']['14N'].q          : eQ            (barn = 10^-28 m^2)\n"
        message += "\t isodata['N']['15N'].conc       : concentration (%/100)\n\n"

        for k,v in self.isodata.items():
            message += f"\t {k} : {v} \n"

        return message


if __name__ == '__main__':

    outcar_spHam = extSpinHamilTensors("./src/OUTCAR_flat",print_io=False)
    #print(outcar_spHam.incar)
    #print(outcar_spHam.outcar)
    print(outcar_spHam)    

