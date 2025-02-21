from spHamilReader.base.pyform import *
from spHamilReader.base.textcolor import C
import os
import sys
import numpy as np
import datetime as dt
import re
#import inspect
#print(f"\n   Code directory : {os.path.dirname(os.path.abspath(__file__))}")
#path=os.path.dirname(os.path.abspath(__file__))
#sys.path.append(path)
#print(sys.path)
'''_________________________________________________________________'''
# Code Information
code_info = \
{
    'summary' : "make AQ tensor file for CCE calculation",
	'mkdate'  : f"{time.ctime(os.path.getctime(__file__))}", 
    'Author'  : "Huijin Park",

    'note1'   : "For VASP, you will find error \n\
    \t\t    if your INCAR file include NGYROMAG\n\
    \t\t    This case need to be considered",

    'note2'   : "POSCAR files should include defect information.\n\
    \t\t    Also, Atomic name need to be consistant",

    'update1' : f"(2023.02.15) Added key : {C.Y}--addhfcore{C.Z} \n\
    \t\t    In the case of VASP, \n\
    \t\t    We should add the core electronic contributions to fermi contact term \n\
    \t\t    (A_tot = A_1c + A_tot)",

    'update2' : f"(2023.03.21) Added key : {C.Y}--fidentify{C.Z} \n\
    \t\t    The first line of the A file include VASP/QUANTUMESSPRESSO.",

    'update3' : f"(2023.03.21) Added key : {C.Y}--version{C.Z} \n\
    \t\t    You can choose version : v1 or v2 \n\
    \t\t    Previous version : Afile of VASP and QE was the same (The central spin factor was treated) \n\
    \t\t    New version : Afile of VASP/QE would be the same to the OUTCAR/hyperfine.out result",

    'update4' : f"(2023.03.21) Added key : {C.Y}--rangeformat{C.Z} \n\
    \t\t    You can convert the format of cell' range from relativepos coordinate to vertex or xyz coordinate",

    'update5' : f"(2023.05.05) Added key : {C.Y}--ccein{C.Z} \n\
    \t\t    You can get interaction tensors of &Exspin in cce.in file. \n\
    \t\t    {C.Y}--JT{C.Z} and {C.Y}--addspins{C.Z} option should be used too.",

    'update6' : f"(2023.05.06) Added key : {C.Y}--cceindiag{C.Z} \n\
    \t\t    If you want to check diagonalized value, use this option. \n\
    \t\t    {C.Y}--JT{C.Z} and {C.Y}--addspins{C.Z} option should be used too.",

    'examples' : "One vacancy defect  : See following example - example01(QE) example02(VASP)\n\
    \t\t    One on-site defect  : See following example - example03(VASP), example05(VASP)\n\
    \t\t    Vacancy cluster     : I have never tried. - try and let me know the issues\n\
    \t\t    Anti-site defect    : I have never tried. - try and let me know the issues\n\
    \t\t    Interstitial defect : I have never tried. - try and let me know the issues"

}
cal_form = cal_form(code_info['summary'],code_info['mkdate'])
################################
if len(sys.argv) == 2:
    if sys.argv[1] == "-h":
        cal_form._start(further_code_info=code_info)
################################

'''_________________________________________________________________'''

# Input parameter Information
argv_info = \
{
'[1]' : '(s)tool           = (s)VASP/QE',
'[2]' : '(s)fposcar        = (s)poscar file',
'[3]' : '(s)dfname         = (s)defect name',
'[4]' : '(s,key)A          = (s)OUTCAR / hyperfine.out , outfilename : Afile',
'[5]' : '(s,key)Q          = (s)OUTCAR / efg.out , outfilename : Qfile',
'[6]' : f'(s,key)D         = (s)OUTCAR / None \n\
(If you only use {C.Y}--D{C.Z} option, put random value in [2]second, [3]third argument)',
'[7]' : '(s,key)dfpos      = (list)defect position x y z (fractional)',
'[8]' : '(s,key)ignoreatm  = (s)ignoring atom (available only for VASP)',
'[9]' : '(s,key)fidentify  = (bool)Add first line or not - VASP orQE [ Opt. True{default} / False ] ',
'[10]': '(s,key)addhfcore  = (bool)Add HF core (A1c) - VASP only [ Opt. True{default} / False ] ',
'[11]': '(s,key)version    = (s)Choose version, default is latest version [ Opt. v1 / v2{default} ] ',
'[12]': '(s,key)rangeformat= (s)The format of defective cell range [ Opt. coord / vertex{default} ]',
'[13]': '(s,key)ccein      = (s)put outfile name - get \"Exspin\" data in cce.in file, electron spin bath case only',
'[14]': f'(i,key)JT         = (int)(IF {C.Y}--ccein(diag){C.Z} is used) Jahn-teller axis ',
'[15]': f'(s,key)addspins   = (s)(IF {C.Y}--ccein(diag){C.Z} is used) sub spin names',
'[16]': '(s,key)cceindiag  = (s)put outfile name - get diagonalized tensors',
}
argv_add_meg = "fposcar, fincar, foutcar : vasp files\
\n** poscar : atoms in vicinity of defect need to be converted in the right atomic name\
\n** poscar : defect atom information is needed"

var = cal_form.is_argv_error(argv_info,argv_add_meg)

'''__________________________________________________________________'''
'''                                                                  '''
'''                      ARGUMENT FUNTION                            '''
'''__________________________________________________________________'''

########################################################################
# Argurment
########################################################################
def argcontrol(var):
    if var.dfpos == []: var.dfpos=None
    else: 
        p = re.compile('[0-9\.]+')
        var.dfpos = list(map(float,p.findall(var.dfpos)))
        print("defect position : ", var.dfpos)
    if var.A==[]: var.A=False
    if var.Q==[]: var.Q=False
    if var.D==[]: var.D=False
    if var.ignoreatm==[]: var.ignoreatm=None
    if var.fidentify=='True' : var.fidentify=True
    if var.fidentify=='False': var.fidentify=False
    if var.addhfcore=='True' : var.addhfcore=True
    if var.addhfcore=='False': var.addhfcore=False
    if var.ccein==[] : var.ccein=False
    if var.cceindiag==[] : var.cceindiag=False
    return var

var = argcontrol(var)
'''__________________________________________________________________'''
'''                                                                  '''
'''                         MAIN PART                                 '''
'''__________________________________________________________________'''

########################################################################
# Funtional Argurment
########################################################################

args = {}

# Quantum Espresso
if var.tool == "QE":

    import spHamilReader.qe as mk

    # "mk" function argument controler
    # Here, the positional argument doesn't need to include "outfile"
    # but, during extrating A,Q,D tensor, the outfile is needed 
    args = \
    {\
    "main" : [var.fposcar,var.dfname],
    "A"    : [var.A ,"Afile"],
    "Q"    : [var.Q ,"Qfile"],
    "D"    : [var.D ,"Dfile"]
    }

# VASP
elif var.tool == "VASP":

    import spHamilReader.vasp as mk

    # "mk" function argument controler
    # Here, the positional argument include "outcar"
    generaloutcar = None

    if var.A : generaloutcar = var.A
    if var.Q : generaloutcar = var.Q
    if var.D : generaloutcar = var.D

    args = \
    {\
    "main" : [generaloutcar, var.fposcar,var.dfname],
    "A"    : ["Afile"],
    "Q"    : ["Qfile"],
    "D"    : ["Dfile"]
    }


else:
    cal_form._end(meg_error=f"Error, tool : '{var.tool}' is wrong value.")

########################################################################
# Main part : Obtain A,Q,D files
########################################################################
do =  mk.do_mkAQDdata(*args["main"],dfposition_frac=None,_ignoreatm=var.ignoreatm,\
                     _fidentify=var.fidentify, _addhfcore=var.addhfcore,_version=var.version,\
                     _rangeformat=var.rangeformat)

if var.A : do.mkAdata(*args["A"])
if var.Q : do.mkQdata(*args["Q"])
#if var.D : do.mkDdata(*args["D"])
if not (var.A or var.Q or var.D) : cal_form._end(meg_error=f"Error, key A,Q,D are not given")

print(do.__repr__([var.A,var.Q,var.D]))

########################################################################
# Postprocess : Obatin cce.in file data 
########################################################################
if var.ccein:

    import spHamilReader.vasp.mkCCEin as mkC

    Adata=None
    Qdata=None
    Ddata=None
    if var.tool=="QE": cal_form._end(meg_error=f"ccein option is available only for VASP")
    if var.JT==[] : cal_form._end(meg_error=f"Error, key JT are not given")
    if var.addspins ==[] : cal_form._end(meg_error=f"Error, key addspins are not given")
    if var.A: Adata = do.getAdata()
    if var.Q: Qdata = do.getQdata()
    if var.D: Ddata = do.getDdata()

    # get data return:data
    data = mkC.dataCCEin()
    data = mkC.getCCEinData(data,do.relativepos['r_1i']\
                 ,var.addspins,var.JT\
                 ,var.A,var.Q,var.D\
                 ,Adata,Qdata,Ddata)

    # read the cce.in if existing
    rdata = mkC.dataCCEin() 
    rdata = mkC.readPreviousFile(rdata,var.ccein)

    # write the data in cce.in
    megs = mkC.megCCEin()
    mkC.writeNewFile(var.ccein,rdata,data,megs,diag=False)

########################################################################
# Postprocess : Obatin ccediag.in file data 
########################################################################
if var.cceindiag:

    import spHamilReader.vasp.mkCCEin as mkC

    Adata=None
    Qdata=None
    Ddata=None
    if var.tool=="QE": cal_form._end(meg_error=f"ccein option is available only for VASP")
    if var.JT==[] : cal_form._end(meg_error=f"Error, key JT are not given")
    if var.addspins ==[] : cal_form._end(meg_error=f"Error, key addspins are not given")
    if var.A: Adata = do.getAdata()
    if var.Q: Qdata = do.getQdata()
    if var.D: Ddata = do.getDdata()

    # get data return:dataDiag
    dataDiag = mkC.dataCCEinDiag()
    dataDiag = mkC.getCCEinDiagData(dataDiag,do.relativepos['r_1i']\
                 ,var.addspins,var.JT\
                 ,var.A,var.Q,var.D\
                 ,Adata,Qdata,Ddata)

    # read the cce.in if existing
    rdataDiag = mkC.dataCCEinDiag()
    rdataDiag = mkC.readPreviousFile(rdataDiag,var.cceindiag)

    # write the data in cce.in
    megsDiag = mkC.megCCEinDiag()
    mkC.writeNewFile(var.cceindiag,rdataDiag,dataDiag,megsDiag,diag=True)

'''__________________________________________________________________'''

cal_form._end()
