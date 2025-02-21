#!/usr/bin/env python
import numpy as np
from spHamilReader.spin_database.constants import * 
from spHamilReader.base.cnvt_coord import frac2cart 
from spHamilReader.base.Calgeometry import distance
from spHamilReader.v1_qe.extract_spinHamils_qe import * 
###############################################################
#
#        FILE FORMAT
#
###############################################################

##########################
# comm
##########################

def divline():
    return  "{:->170}\n".format("")

##########################
# A tensor file
##########################

def Afilehead():
    return   "{:>10}\t{:>10}\t{:>10}\t{:>10}\t\
{:>10}\t{:>10}\t{:>10}\t{:>10}\t\
{:>10}\t{:>10}\t{:>10}\t\
{:>10}\t{:>10}\t{:>10}\n"\
.format("AtomName","r_1i_X","r_1i_Y","r_1i_Z",\
        "iso[MHz]", "Axx", "Axy", "Axz",\
                    "Ayx", "Ayy", "Ayz",\
                      "Azx", "Azy", "Azz")                 

def Afileinfo():
    atm = "{:5}"
    gfactorline = "g-factor___: {:>10f} [No-unit]\n"
    etcline     = "etc____:\t{:10}\t{:10}\t{:10}\t{:10.5f}\n"\
                              .format("---","---","---",0.00)
    minline = "MinDif[A]:\t{:10.5f}\t{:10.5f}\t{:10.5f}\n"
    maxline = "MaxDif[A]:\t{:10.5f}\t{:10.5f}\t{:10.5f}\n"

    Dict = {'gline'     : atm+gfactorline,
            'etcline'     : atm+etcline,
            'minline'    : minline,
            'maxline'    : maxline}

    return Dict

def Atensorform():
    return "{:10}\t{:10.5f}\t{:10.5f}\t{:10.5f}\t\
{:10.5f}\t{:10.5f}\t{:10.5f}\t{:10.5f}\t\
{:10.5f}\t{:10.5f}\t{:10.5f}\t\
{:10.5f}\t{:10.5f}\t{:10.5f}\n"


##########################
# Q tensor file
##########################

def Qfilehead():
    return   "{:>10}\t{:>10}\t{:>10}\t{:>10}\t\
{:>10}\t{:>10}\t{:>10}\t\
{:>10}\t{:>10}\t{:>10}\t\
{:>10}\t{:>10}\t{:>10}\t\
[Hatree/Bohr_radius^2]\n"\
.format("AtomName","r_1i_X","r_1i_Y","r_1i_Z"\
                      ,"EFGxx","EFGxy","EFGxz"\
                      ,"EFGyx","EFGyy","EFGyz"\
                      ,"EFGzx","EFGzy","EFGzz")    

def Qfileinfo():
    atm = "{:5}"
    eQline         = "eQ___:\t{:>10f}\t[Q/milibarn]*(10^-1)\n"
    etcline     = "etc____:\t{:10}\t{:10}\t{:10}\t"\
                              .format("---","---","---")
    etclinetail= "{:10.5f}\t{:10.5f}\t{:10.5f}\t\
{:10.5f}\t{:10.5f}\t{:10.5f}\t\
{:10.5f}\t{:10.5f}\t{:10.5f}\n"
    minline = "MinDif[A]:\t{:10.5f}\t{:10.5f}\t{:10.5f}\n"
    maxline = "MaxDif[A]:\t{:10.5f}\t{:10.5f}\t{:10.5f}\n"

    Dict = {'eQline'     : atm+eQline,
            'etcline'     : atm+etcline+etclinetail,
            'minline'    : minline,
            'maxline'    : maxline}

    return Dict

def Qtensorform():
    return "{:10}\t{:10.5f}\t{:10.5f}\t{:10.5f}\t\
{:10.5f}\t{:10.5f}\t{:10.5f}\t\
{:10.5f}\t{:10.5f}\t{:10.5f}\t\
{:10.5f}\t{:10.5f}\t{:10.5f}\n"

def addelement(_list,elem):
    if elem not in _list:
        _list.append(elem)
    return _list

def addlist(_list,sublist):
    if sublist not in _list:
        _list.append(sublist)
    return _list

###############################################################
#
#        Check poscar if the file include central defect
#
###############################################################
def doesItHaveDefect(dposcar,dfctname):
    if dfctname in dposcar["atom_data"]:
        return  dposcar["atom_data"].index(dfctname)
    else:
        return None

def mkDefectKeyInPOSCAR(dposcar,dfctname,defectposition_frac=None):

    """ To simplify, make a key for only central defect in poscar data"""

    idx = doesItHaveDefect(dposcar,dfctname)

    if idx != None: 
        # remove defect atomic position in unit cell information
        defectindex = sum(dposcar['atom_data_num'][:idx])
        defectposition = dposcar['unit_cell_information'][defectindex]
        del dposcar['unit_cell_information'][defectindex]
        del dposcar['atom_data_num'][idx]
        del dposcar['atom_data'][idx]

        # add defect atomic position as a new key
        dposcar['central_defect'] = defectposition    
    else:
        if defectposition==None:
            sys.exit("Error!!! \
                        you should add defect position in argument of code \
                        or add defect position in file")
        else:
            dposcar['central_defect'] = [dfctname] + defectposition    

    return dposcar    

###############################################################
#
#         Relative Distance 
#
###############################################################

def relativeDistance(latt_constant,unit_cell_parameter,frac_unitCellInfo,frac_centralDefect):
    
    relativepos = { 'r_1i': [],
                    'MinDif[A]': [], 
                    'MaxDif[A]': [],
                    'farthest_atoms_idx' : [],
                    'farthest_atoms_info' : []} 

    cellparam = np.asarray(unit_cell_parameter)*latt_constant
    cart_centralDefect = frac2cart(cellparam,\
                            np.array(frac_centralDefect[1:4]))

    r_1i_dist = []

    atm_name=frac_unitCellInfo[0][0]

    for i,frac_atm in enumerate(frac_unitCellInfo):

        cart_atm = frac2cart(cellparam,\
                        np.array(frac_atm[1:4]))

        # Find farthest atom to etc line ##########################
        r_1i_dist.append([i,distance(cart_atm,cart_centralDefect)])
        if atm_name != frac_atm[0] or i == len(frac_unitCellInfo)-1:
            if i != len(frac_centralDefect)-1:
                del r_1i_dist[-1]

            zip_r_1i_dist = list(zip(*r_1i_dist))
            idx = zip_r_1i_dist[1].index(max(zip_r_1i_dist[1]))
            farthest_atom_idx  = r_1i_dist[idx][0]
            farthest_atom_info = relativepos['r_1i'][farthest_atom_idx]
            relativepos['farthest_atoms_idx'].append(farthest_atom_idx)
            relativepos['farthest_atoms_info'].append(farthest_atom_info)
            r_1i_dist = []
            r_1i_dist.append([i,distance(cart_atm,cart_centralDefect)])
        atm_name= frac_atm[0]
        ###########################################################

        r_1i = cart_atm - cart_centralDefect

        relativepos['r_1i'].append( [frac_atm[0]] + list(r_1i) )

        if frac_atm[0] == frac_centralDefect[0]:# atomic name..
            sys.exit("Error!!! \
                        There is central defect in unit cell information")

    ri_zip = list(zip(*relativepos['r_1i']))

    relativepos['MinDif[A]'] = [min(ri_zip[1]), min(ri_zip[2]), min(ri_zip[3])]
    relativepos['MaxDif[A]'] = [max(ri_zip[1]), max(ri_zip[2]), max(ri_zip[3])]

    return relativepos 

###############################################################
#
#        MAKE A TENSOR FILE
#
###############################################################

def Afile(relativepos,doutcar_Adata,doutcar_isodata,fout):

    """
        Args:

            relativepos(dict)
             - Keys 
                 'r_1i'(2d list)     : relative atomic position 
                                     between atom and central defect 
                                     e.g. [Xi,Yi,Zi] - [X_df,Y_df,Z_df]
                'MinDif[A]'(list) : minimum r_1i
                'MaxDif[A]'(list) : maximum r_1i

            doutcar_Adata(dict) 
             - Keys 
                'Aiso'            :[]        # ion# x Fermi contact(MHz) 
                'Adip'            :[]        # ion# x [Axx, Axy, Axz, Ayx .., Azz] (MHz)
                'Adiag'            :[]        # ion# x [Axx(D), Ayy(D), Azz(D)] (MHz) for total A
                  'giveng'        :[]        # ion# x given g factor (μ/μN) in gipaw.in

            doutcar_isodata(dict)
             - Keys
                key(str) : atom name (ex. B)
                value(myspindict) : isotope information for each atom (ex. Isodata['B'] =  {'10B': (.. , .., ) })
                ==> ex. : Isodata['B'] =  {'10B': (3.0, 2.8747, 0.08450, 0.1990) }
                isodata['B']['10B'].s          : spin         (no unit)
                isodata['B']['11B'].gyro      : gyro..     (rad/ms/G)
                isodata['N']['14N'].q          : eQ         (barn = 10^-28 m^2)
                isodata['N']['15N'].conc      : concent.. (%/100)

            fout(str)
    """

    dgfactor=[]
    dmindiff=[]
    dmaxdiff=[]
    dreldist=[]
    dAtensor=[]
    
    # g factor 
    # NOTE     : nuclear magnetic moment(\mu) = g * nuclear magneton(\mu_N) * I / h_bar
    #        : g = ( h_bar * \mu ) / ( I * \mu_N )
    #        : gyro magnetic ratio(\gamma_n) = g * nuclear magneton(\mu_N) / h_bar
    #        : g = ( h_bar * \gamma_nn ) / \mu_N

    dgfactor = []
    for atm, isotmp in doutcar_isodata.items():
        for iso, dat in isotmp.items(): 
            g = dat.gyro * 1e7 * HBAR_SI / NUCLEAR_MAGNETON
            dgfactor = addlist(dgfactor,[iso,g])

    Dict = Afileinfo()
    # Make A tensor file
    with open(fout, 'w') as f:
        # HEAD
        f.write("QUANTUMESSPRESSO\n")
        f.write(divline())
        f.write(Afilehead())
        # ISOTOPE DATA
        f.write(divline())
        for iso_g in dgfactor:
            f.write(Dict['gline'].format(*iso_g))
        for iso_g in dgfactor:
            f.write(Dict['etcline'].format(iso_g[0]))
        f.write(Dict['minline'].format(*relativepos['MinDif[A]']))
        f.write(Dict['maxline'].format(*relativepos['MaxDif[A]']))
        # MAIN DATA
        f.write(divline())
        for i,atm in enumerate(relativepos['r_1i']):
            corrected_Aiso = np.asarray(doutcar_Adata['Aiso'][i])/np.asarray(doutcar_Adata['giveng'][i])
            corrected_Adip = np.asarray(doutcar_Adata['Adip'][i])/np.asarray(doutcar_Adata['giveng'][i])
            f.write(Atensorform().format(*atm,\
                                        corrected_Aiso,\
                                        *corrected_Adip))


def Qfile(relativepos,doutcar_Qdata,doutcar_isodata,fout,ignoreatm=None):

    # Already : Hatree/(Bohr_radius^2) 

    # eQ [unit : barn = 10^-28 m^2 ]

    deQ = []
    for atm, isotmp in doutcar_isodata.items():
        for iso, dat in isotmp.items(): 
            # barn to "milibarn(10^-31 m^2) * 10^1" (10^-30 m^2)
            eQ = dat.q * pow(10,2) # * 10^-2 (10^-28 m^2) = (10^-30 m^2)
            deQ = addlist(deQ, [iso,eQ])

    # Farthest atoms EFG tensor
    detcEFG = []
    i=0
    for atm, isotmp in doutcar_isodata.items():
        if atm != ignoreatm :
            idx = relativepos['farthest_atoms_idx'][i]
            EFG = np.array(doutcar_Qdata['EFG'][idx]) 
            i += 1
            for iso, dat in isotmp.items():
                detcEFG = addlist(detcEFG,[iso, EFG])

    Dict = Qfileinfo() 
    # Make Q tensor file
    with open(fout, 'w') as f:
        # HEAD
        f.write(divline())
        f.write(Qfilehead())
        # ISOTOPE DATA
        f.write(divline())
        for iso_eQ in deQ:
            f.write(Dict['eQline'].format(*iso_eQ))
        for iso_etcEFC in detcEFG:
            f.write(Dict['etcline'].format(iso_etcEFC[0],*iso_etcEFC[1]))

        f.write(Dict['minline'].format(*relativepos['MinDif[A]']))
        f.write(Dict['maxline'].format(*relativepos['MaxDif[A]']))

        # MAIN DATA
        f.write(divline())
        for i,atm in enumerate(relativepos['r_1i']):
            EFG = np.array(doutcar_Qdata['EFG'][i])
            f.write(Qtensorform().format(*atm,*EFG))

class do_mkAQDdata:

    def __init__(self,fposcar,dfname,dfposition_frac=None,_ignoreatm=None):
        self.fposcar = fposcar
        self.dfname = dfname
        self._ignoreatm = _ignoreatm
        self.dfposition_frac = dfposition_frac

        print(self.startmessage())
        #print(self.fposcar)
        #print(self.dfname)
        #print(self.dfposition_frac)
        # df position : fractional ..

        ###########################################################
    
        # Calculate relative position 
    
        from ioPoscar import read_POSCAR
        from dictproc import print_dict
        self.poscar = read_POSCAR(fposcar)
        self.poscarmod = mkDefectKeyInPOSCAR(self.poscar,dfname,defectposition_frac=dfposition_frac)
        #if you want to make central defect
        #    or the defect information is not in poscar data
        #     , you should add defect position 
        #To obtain relative position, defect position is needed
        self.relativepos  = relativeDistance(self.poscarmod["latt_constant"],\
                                             self.poscarmod["unit_cell_parameter"],\
                                             self.poscarmod["unit_cell_information"],\
                                             self.poscarmod['central_defect'])
    
        ###########################################################

        # Read outfiles informations (A,Q,D tensor,,, isotopes... etc)
        self.ext = extSpinHamilTensors(self.poscarmod)
        self.isodata = (self.ext).isodata
        ###########################################################

    def mkAdata(self,fhfout,fout_Afile):
        print("\n")
        self.fhfout = fhfout
        self.fout_Afile = fout_Afile
        (self.ext).read_Atensor(self.fhfout)
        self.Adata = (self.ext).Adata

        #print("\t",self.fhfout)
        #print("\t",self.fout_Afile)

        # Make A tensor file
        Afile(self.relativepos, self.Adata, self.isodata,self.fout_Afile)
        ###########################################################

    def mkQdata(self,fefgout,fout_Qfile):
        print("\n")
        self.fefgout = fefgout    
        self.fout_Qfile = fout_Qfile
        (self.ext).read_Qtensor(self.fefgout)
        self.Qdata = self.ext.Qdata

        #print("\t",self.fefgout)
        #print("\t",self.fout_Qfile)

        # Make A tensor file
        Qfile(self.relativepos, self.Qdata, self.isodata,self.fout_Qfile,ignoreatm=self._ignoreatm)
        ###########################################################

    def startmessage(self):
        message="\n\t=========================================\n"
        message+="\t START TO WRITE SPIN TENSOR FILES \n"
        message+="\t========================================="
        return message

    def __repr__(self):
        message = ""
        message = "\n"
        message += "\t" + "-"*40 + "\n"
        message += f"\t{type(self).__name__} : \n"
        message += "\t" + "-"*40 + "\n"

        message+="\n\t=========================================\n"
        message+="\t READ OUTCAR (KEY : outcar_spHam) \n"
        message+="\t=========================================\n"
        message+=self.ext.__repr__()
        message+="\n\t=========================================\n"
        message+="\t READ POSCAR (KEY : poscar or poscarmod) \n"
        message+="\t=========================================\n"
        message+=self.prt_poscar()
        message+="\n\t=========================================\n"
        message+="\t          MAKE A,Q FILE                     \n"
        message+="\t=========================================\n"
        message+=self.prt_relativepos()
        message+=f"\n"

        try:
            self.Adata['Adip']
            message+=f"\tAfile : {fout_Afile}\n"
        except:
            message+=f"\tAfile : No read A\n"
            
        try:
            self.Qdata['EFG'] 
            message+=f"\tQfile : {fout_Qfile}\n"
        except:
            message+=f"\tQfile : No read Q\n"

        return message

    def prt_poscar(self):
        message = "" 
        message += f"\n\t POSCAR from {self.fposcar} (KEY : poscar)\n"
        message += "\t" + "-"*30 + "\n"
        message += f"\t *atom_data : {self.poscar['atom_data']}\n" 
        message += f"\t *atom_data_num : {self.poscar['atom_data_num']}\n"
        message += f"\t *unit_cell_information :  \n\n"
        message += printform(3,["frac_x","frac_y","frac_z"],\
                                [v[1:] for v in self.poscar["unit_cell_information"]])
        return message

    def prt_relativepos(self):
        message = "" 
        message += f"\n\t Relative positions from central defect (KEY : relativepos)\n"
        message += "\t" + "-"*30 + "\n"
        message += f"\t *MinDif[A] : {self.relativepos['MinDif[A]']}\n"
        message += f"\t *MaxDif[A] : {self.relativepos['MaxDif[A]']}\n"
        message += f"\t *farthest_atoms_idx : {self.relativepos['farthest_atoms_idx']}\n"
        message += f"\t *farthest_atom_info : {self.relativepos['farthest_atoms_info']}\n"
        message += f"\t *r_1i (cart) :  \n\n"
        message += printform(4,["atom","r_1ix","r_1iy","r_1iz"],\
                                [v[:] for v in self.relativepos['r_1i']])
        return message



if __name__ =='__main__':

#    ###########################################################
#    
#    # Read outcar informations (A,Q,D tensor,,, isotopes... etc)
#
#    from extract_spinHamils import extSpinHamilTensors
#    outcar_spHam = extSpinHamilTensors("./src/OUTCAR_flat",print_io=False)
#    #print(outcar_spHam)
#    ###########################################################

    ###########################################################

#    ###########################################################
#    # Example : 
#    ###########################################################
#
#
#    fposcar = "./src/QE_POSCAR_flat"
#    fhfout  = "./src/QE_Hyper_flat.out"
#    fefgout  = "./src/QE_EFG_flat.out"
#
#    dfname  = 'VB'
#    fout_Afile = "./src/QE_ATensor_flat"
#    fout_Qfile = "./src/QE_QTensor_flat"
#    do =  do_mkAQDdata(fposcar,dfname,dfposition_frac=None)
#    do.mkAdata(fhfout,fout_Afile)
#    do.mkQdata(fefgout,fout_Qfile)
#    print(do)

#    ###########################################################
#    # Main : 
#    ###########################################################
#    import sys
#    if len(sys.argv)!=7:
#        print("[1] hyper.out")
#        print("[2] efg.out")
#        print("[3] fposcar")
#        print("[4] dfname")
#        print("[5] fout_Afile")
#        print("[6] fout_Qfile")
#        sys.exit('Check the arguments')
#
#    fhfout  = sys.argv[1]
#    fefgout = sys.argv[2]
#    fposcar = sys.argv[3]
#    dfname  = sys.argv[4]
#    fout_Afile = sys.argv[5]
#    fout_Qfile = sys.argv[6]
#
#    do =  do_mkAQDdata(fposcar,dfname,dfposition_frac=None)
#    do.mkAdata(fhfout,fout_Afile)
#    do.mkQdata(fefgout,fout_Qfile)
#    print(do)

    ###########################################################
    # Only A 
    ###########################################################


    fposcar = "./src/QE_POSCAR_flat"
    fhfout  = "./src/QE_Hyper_flat.out"

    dfname  = 'VB'
    fout_Afile = "./a"#"./src/QE_ATensor_flat"
    do =  do_mkAQDdata(fposcar,dfname,dfposition_frac=None)
    do.mkAdata(fhfout,fout_Afile)
    print(do)


