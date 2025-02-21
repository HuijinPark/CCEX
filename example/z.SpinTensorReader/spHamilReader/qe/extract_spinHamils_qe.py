import numpy as np
from spHamilReader.base.rwfiles import readfile
from spHamilReader.qe.read_qe import *
from spHamilReader.base.printatm import *
from spHamilReader.base.ioPoscar import read_POSCAR
from spHamilReader.base.dictproc import print_dict
from spHamilReader.spin_database.loadisotopes import isotopeData

def readtensor_3x3(linedata1,linedata2,linedata3):
    tensor1 = list(map(float,linedata1))
    tensor2 = list(map(float,linedata2))
    tensor3 = list(map(float,linedata3))
    return    [tensor1[0],tensor1[1],tensor1[2],\
             tensor2[0],tensor2[1],tensor2[2],\
             tensor3[0],tensor3[1],tensor3[2]]

def readtensor_diag(linedata):
    tensor = list(map(float,linedata))
    return    [tensor[0],tensor[1],tensor[2]]

def readtensor_elements(linedata,elem):
    tensor = list(map(float,linedata))
    return    [tensor[i] for i in range(elem)]

def readtensor_Aiso(linedata):
    return float(linedata[3])    

def readtensor_given_gfactor(linedata):
    return float(linedata[0])    

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

    keys={'EFG'            :[],    # ion# x [Vxx, Vxy, Vxz, Vyx .., Vzz] (V/A^2)
          'EFG_diag'    :[],    # ion# x [Vxx(D), Vyy(D), Vzz(D)] (V/A^2) 
          'Eigenvectors':[],    # ion# x [vxx, vxy, vxz, vyx .., vzz]
          'Cq'            :[],    # ion# x Cq (MHz)
          'eta'            :[],    # ion# x eta 
          'q'            :[],    # ion# x Q (mb)
          'Cq,eta,q'    :[]}
    
    ###############################################
    read=False
    line=startingline
    while(line<len(dat)):
        ###############################################
        # find EFG tensor
        if checkLine_all(dat[line],"----- total EFG (symmetrized) -----"):
            readtensor=readtensor_3x3
            key='EFG'
            read=True
        ## find EFG tensor diag
        #if checkLine_all(dat[line]," Electric field gradients after diagonalization (V/A^2)"):
        #    readtensor=readtensor_diag
        #    key='EFG_diag'
        #    read=True
        ## find Cq, eta,q
        #if checkLine_all(dat[line],"NMR quadrupolar parameters"):
        #    readtensor=readtensor_diag
        #    key='Cq,eta,q'
        #    read=True
        ###############################################
        if read == True:
            nion=0
            while(nion<atomicdata["total_atom_num"]):
                try:
                    if key == 'EFG':
                        keys[key].append(readtensor(dat[line][2:],dat[line+1][2:],dat[line+2][2:]))
                        line+=3
                    nion +=1
                except: pass;
                line+=1
            read=False
        line+=1
        ###############################################
        #break while when you got all materials
        #if keys['Cq,eta,q'] != []:
        #    break;
    ###############################################

    #tmp = np.array(keys['Cq,eta,q'])
    #keys['Cq']  = tmp[:,0]
    #keys['eta'] = tmp[:,1]
    #keys['q']   = tmp[:,2]
    #del keys['Cq,eta,q']

    if print_io == True:
        for k,v in keys.items():    
            print("-"*40)
            if k=='EFG' or k == 'EFG_diag':
                print(k,'(Hartree/bohrradius^2)')
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

    return keys

def extAdata(atomicdata,dat,print_io=False,startingline=0):

    keys={'Aiso'    :[],    # ion# x Fermi contact(MHz) 
          'Adip'    :[],    # ion# x [Axx, Axy, Axz, Ayx .., Azz] (MHz)
          'Adiag'   :[],    # ion# x [Axx(D), Ayy(D), Azz(D)] (MHz) for total A
          'giveng'  :[]}    # ion# x given g factor (μ/μN) in gipaw.in
    ###############################################
    inputtedgfactor=[]
    read=False
    line=startingline
    while(line<len(dat)):
        ###############################################
        # find g factor given as input value
        if checkLine_all(dat[line],\
            "NUCLEAR G-TENSORS FROM INPUT:"):
            line +=1
            readtensor=readtensor_given_gfactor
            key='giveng'
            read=True
        # find Fermi contact 
        if checkLine_all(dat[line],\
            "----- Fermi contact in MHz ----- "):
            line +=1
            readtensor=readtensor_Aiso
            key='Aiso'
            read=True
        # find A dipolar term 
        if checkLine_all(dat[line],\
            "----- total dipolar (symmetrized) -----"):
            readtensor=readtensor_3x3
            key='Adip'
            read=True
        # For qe data, there are no diagonalize data for full tensor
        ###############################################
        if read == True:
            nion=0
            while(nion<atomicdata["total_atom_num"]):
                try:
                    if key == 'Aiso':
                        keys[key].append(readtensor(dat[line][2:]))
                    elif key == 'Adip':
                        keys[key].append(readtensor(dat[line][2:],dat[line+1][2:],dat[line+2][2:]))
                        line+=3
                    elif key == 'giveng':
                        keys[key].append(readtensor(dat[line][2:]))
                    nion +=1
                except: pass;
                line+=1
            read=False
        line+=1
        ###############################################
        #break while when you got all materials
        #if keys['Adiag'] != []:
        #    break;
    ###############################################

    if print_io == True:
        for k,v in keys.items():    
            print("-"*40)
            if k != 'giveng':
                print(k,"(MHz)")
            else:
                print(k,"(μ/μN)")
            print("-"*40)
            for j in range(len(v)):
                print(j, v[j])
            print("-"*40)

    return keys

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
                'EFG'            :[]        # ion# x [Vxx, Vxy, Vxz, Vyx .., Vzz] (V/A^2)
                'EFG_diag'        :[]        # ion# x [Vxx(D), Vyy(D), Vzz(D)] (V/A^2) 
                'Cq'            :np.arr    # ion# x Cq (MHz)
                'eta'            :np.arr    # ion# x eta 
                'q'                :np.arr    # ion# x Q (mb)
                ## None ! 'Eigenvectors' :[] # ion# x [vxx, vxy, vxz, vyx .., vzz] 

            For A tensor : Adata
                'Aiso'            :[]        # ion# x Fermi contact(MHz) 
                'Adip'            :[]        # ion# x [Axx, Axy, Axz, Ayx .., Azz] (MHz)
                'Adiag'            :[]        # ion# x [Axx(D), Ayy(D), Azz(D)] (MHz) for total A
                  'giveng'        :[]        # ion# x given g factor (μ/μN) in gipaw.in

            For D tensor : Ddata
                'ZFS'            :np.arr    # np.array([Dxx, Dxy, Dxz, Dyx .., Dzz]) (MHz)
                'D'                :float    # 3/2*Dzz (MHz)
                'E'                :float    # (Dxx - Dyy)/2 (MHz)
                'D_diag'        :[]        # [Dxx(D),Dyy(D),Dzz(D)] (MHz) 
                'D_eigvec'        :[]        # [[vec_x, vec_y, vec_z],[],[]]

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


    def __init__(self,fposcar):

        if type(fposcar)==str:
            self.atomicdata = read_POSCAR(fposcar)
        else:
            self.atomicdata = fposcar 

        self.atomicdata['total_atom_num'] = len(self.atomicdata['unit_cell_information'])
        #print_dict(self.atomicdata)
        self.Qdata={}
        self.Adata={}
        self.Ddata={}
        self.isodata={}

        c=0
        #print(self.atomicdata['atom_data'])
        #print(isotopeData[self.atomicdata['atom_data'][0]])
        for v in self.atomicdata['atom_data']:
            #print("atm:",v)
            #print("isotopeData",isotopeData(v))
            try:
                self.isodata[v] = isotopeData(v) 
                c+=1
            except:
                for i in range(c):
                    #print(i,self.atomicdata['atom_data'][i])
                    if self.atomicdata['atom_data'][i] in v:
                        self.isodata[v] = self.isodata[self.atomicdata['atom_data'][i]]
                        break;
            #print("iso",v,self.isodata[v])

    def read_Atensor(self,fhfout,print_io=False):

        print(f"\tRead hyperfine tensor from {fhfout} ...")
        self.dhf = readfile(fhfout)
        self.Adata = extAdata(self.atomicdata,self.dhf\
                            ,print_io=print_io)

    def read_Qtensor(self,fefgout,print_io=False):

        print(f"\tRead efg tensor from {fefgout} ...")
        self.defg = readfile(fefgout)
        self.Qdata= extQdata(self.atomicdata,self.defg\
                            ,print_io=print_io)

        
    def __repr__(self):
        message = "\n"
        message += "\t"+ "-"*40 + "\n"
        message += "\t"+f"{type(self).__name__} : \n"
        message += "\t"+"-"*40 + "\n"

        message += self.prt_atomicdata()

        try:
            self.Qdata['EFG']
            message += self.prt_readQ()
        except:
            pass;
        
        try: 
            self.Adata['Adip']
            message += self.prt_readA()
        except:
            pass;

        #if 'T' in self.incar['LDMATRIX']:
        #    message += self.prt_readD()
        
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
        message += "\n\t*EFG : EFG tensor (Hartree/bohrradius^2) \n\n"
        message += printform(9,[ "Vxx","Vxy","Vxz",\
                                 "Vyx","Vyy","Vyz",\
                                 "Vzx","Vzy","Vzz"],EFG)
        ###################################
        #message += "\n*EFG_diag : Diagonalized EFG tensor (V/A^2) \n\n"
        #message += printform(3,["Vxx(D)","Vyy(D)","Vzz(D)"],EFG_diag)
        ###################################
        #message += "\nNMR quadrupolar parameters\n\n"
        #message += "*Cq  : quadrupolar parameter    Cq=e*Q*V_zz/h (MHz) \n"
        #message += "*eta : asymmetry parameters     (V_yy - V_xx)/ V_zz \n"
        #message += "*Q   : nuclear electric quadrupole moment in mb (millibarn = 10^-31 m^2) \n\n"
        #message += printform(3,["Cq","eta","q"],list(zip(Cq,eta,q)))
        message += "\t"+"-"*30 + "\n"
        ###################################
        return message

    def prt_readA(self):
        Aiso = self.Adata['Aiso']
        Adip = self.Adata['Adip']
        Adiag = self.Adata['Adiag']
        giveng = self.Adata['giveng']
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
        message += "\n\t*giveng : given g factor in gipaw.in file  (μ/μN) \n"
        message += "\n\t  ( #NOTE A tensor file will be written\n\
          after dividing A tensor into giveng values ) \n\n"
        message += printform(1,["g_factor"],[[g] for g in giveng])
        ###################################
        #message += "\n*Adiag : Diagonalized total Hyperfine tensor (MHz) \n\n"
        #message += printform(3,[ "Axx(D)","Ayy(D)","Azz(D)"],Adiag)
        message += "\t" + "-"*30 + "\n"
        ###################################
        return message

#    def prt_readD(self):
#        message = "" 
#        message += "\nZero Field splitting tensor (KEY : Ddata)\n"
#        message += "-"*30 + "\n"
#
#        message += "\n*ZFS : ZFS tensor (MHz) \n\n"
#        message += "{:>10.3f}{:>10.3f}{:>10.3f}\n".format(self.Ddata['ZFS'][0,0],self.Ddata['ZFS'][0,1],self.Ddata['ZFS'][0,2])
#        message += "{:>10.3f}{:>10.3f}{:>10.3f}\n".format(self.Ddata['ZFS'][1,0],self.Ddata['ZFS'][1,1],self.Ddata['ZFS'][1,2])
#        message += "{:>10.3f}{:>10.3f}{:>10.3f}\n".format(self.Ddata['ZFS'][2,0],self.Ddata['ZFS'][2,1],self.Ddata['ZFS'][2,2])
#            
#        message += "\nDiagonalized ZFS tensor \n"
#        message += "*D_diag : Diagonalized ZFS tensor (MHz) \n"
#        message += "*D_eigvec : Eigenvectors for diagonalized ZFS tensor \n\n"
#        message += printform(4,[ "Dii(D)", 'eigvec_x', 'eigvec_y', 'eigvec_z'],\
#            [[d]+v for d,v in zip(self.Ddata['D_diag'],self.Ddata['D_eigvec'])])
#
#        message += "\nZFS parameters (MHz) \n"
#        message += "*D : axial component D = 3/2*D_zz (MHz) \n"
#        message += "*E : transverse component E = (D_xx-D_yy)/2 (MHz) \n\n"
#        message += "{:>7}{:>10.3f}".format('D = ',self.Ddata['D']) + "\n"
#        message += "{:>7}{:>10.3f}".format('E = ',self.Ddata['E']) + "\n"
#        message += "-"*30 + "\n"
#        return message

    def prt_atomicdata(self):
        message = "" 
        message += "\n\tAtomic data from POSCAR (KEY : atomicdata)\n"
        message += "\t" + "-"*30 + "\n"

        message += f"\t *atom_data : {self.atomicdata['atom_data']}\n" 
        message += f"\t *total_atom_num : {self.atomicdata['total_atom_num']}\n"
        message += f"\t *unit_cell_information :  \n\n"
        message += printform(3,["frac_x","frac_y","frac_z"],\
                                [v[1:] for i,v in enumerate(self.atomicdata["unit_cell_information"])])
        return message

    def prt_isotopedata(self):
        message = "" 
        message += "\n\tIsotopic data (KEY : isodata)\n"
        message += "\t" + "-"*30 + "\n\n"

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

    ext = extSpinHamilTensors("./src/POSCAR_strn100")
    ext.read_Atensor("./src/hfn100.out",print_io=False)
    ext.read_Qtensor("./src/efg100.out",print_io=False)
    print(ext)
