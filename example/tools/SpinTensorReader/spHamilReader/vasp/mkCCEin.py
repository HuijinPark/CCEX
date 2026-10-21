import re
from spHamilReader.base.rwfiles import readfile
import numpy as np

def getSpinNames(varAddspins):
    # check addspin
    p = re.compile('[0-9a-zA-Z]+')
    return list(p.findall(varAddspins))

def megCCEin():
    # meg
    megs = \
    {
    "eq"  : "eqs = {:>5s} {:>10.5f}\n",
    "rxyz": "rxyz({:d}) = {:>5s} {:>10.5f} {:>10.5f} {:>10.5f}\n",
    "A"   : "A({:d}) = {:>5s} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f}\n",
    "Q"   : "Q({:d}) = {:>5s} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5f}\n",
    "D"   : "ZFS({:d}) = {:>5s} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f}\n"
    }
    return megs

def megCCEinDiag():
    # meg
    megs = \
    {
    "eq"  : "eqs = {:>5s} {:>10.5f}\n",
    "rxyz": "rxyz({:d}) = {:>5s} {:>10.5f} {:>10.5f} {:>10.5f}\n",
    "Adiag": "Adiag({:d}) = {:>5s} {:>10.5f} {:>10.5f} {:>10.5f}\n",
    "Qdiag": "Qdiag({:d}) = {:>5s} {:>10.5f} {:>10.5f} {:>10.5f}\n",
    "Ddiag": "Ddiag({:d}) = {:>5s} {:>15.5f} {:>15.5f} {:>15.5f}\n",
    }
    return megs

def dataCCEin():
    # data 
    data = \
    {
    "eq"  : [], 
    "rxyz": [],
    "A"   : [],
    "Q"   : [],
    "D"   : [],
    }
    return data

def dataCCEinDiag():
    # data 
    data = \
    {
    "eq"  : [], 
    "rxyz": [],
    "Adiag": [],
    "Qdiag": [],
    "Ddiag": [],
    }
    return data

def appendeq(obj,addspin,eq):
    return obj + [[addspin,eq]]
def appendrxyz(obj,jt,addspin,atm):
    return obj + [[jt,addspin,atm[1],atm[2],atm[3]]]
def appendtensor(obj,jt,addspin,tensor):
    return obj + [[jt,addspin,tensor[0,0],tensor[0,1],tensor[0,2],tensor[1,0],tensor[1,1],tensor[1,2],tensor[2,0],tensor[2,1],tensor[2,2]]]
def appendtensordiag(obj,jt,addspin,tensor):
    return obj + [[jt,addspin,tensor[0],tensor[1],tensor[2]]]

def setAtensor(Adata,addspin,i):
    return Adata['dgfactor'][addspin] * (np.eye(3) * Adata['cAiso'][i]) \
       + Adata['dgfactor'][addspin] * Adata['cAdip'][i].reshape(3,3)

def setQtensor(Qdata,i):
    return Qdata['cEFG'][i].reshape(3,3)

def setDtensor(Ddata):
    return Ddata['ZFS'].reshape(3,3)

def getCCEinData(data,rxyzs,addspins,jt,Aio,Qio,Dio,Adata,Qdata,Ddata):
    # check addspin
    addspins=getSpinNames(addspins)
    # get A,Q data
    for i,atm in enumerate(rxyzs):
        for j,addspin in enumerate(addspins):
            val = re.findall('[a-zA-Z]+',addspin) # spin name
            if atm[0] == val[0]:
                data['rxyz'] = appendrxyz(data['rxyz'],jt,addspin,atm)
                if Aio:
                    A = setAtensor(Adata,addspin,i)
                    data['A'] = appendtensor(data['A'],jt,addspin,A)
                if Qio:
                    Q  = setQtensor(Qdata,i) 
                    data['Q'] = appendtensor(data['Q'],jt,addspin,Q)
                    data['eq'] = appendeq(data["eq"],addspin,Qdata['ceQ'][addspin])

    # get D data
    if Dio:
        D  = setDtensor(Ddata) 
        data['D'] = appendtensor(data['D'],jt,'e',D)

    return data

def sortEigXYZ_3x3(eigenvalues,eigenstates):
    idx = np.argsort(np.abs(eigenvalues))
    return eigenvalues[idx]
    ## find the main principal axis
    ## e.g.) principalaxis0 : main principal axis of eigenstates[:,0]
    #principalaxis0 = np.where(np.abs(eigenstates[:,0]) == max(np.abs(eigenstates[:,0])))[0]
    #principalaxis1 = np.where(np.abs(eigenstates[:,1]) == max(np.abs(eigenstates[:,1])))[0]
    #principalaxis2 = np.where(np.abs(eigenstates[:,2]) == max(np.abs(eigenstates[:,2])))[0]
    #principalaxes = np.array([principalaxis0,principalaxis1,principalaxis2]).reshape(3,)

    ## Check if maxidx1,2,3 are one of the value 0(x),1(y),2(z)
    #if np.dot(principalaxes,principalaxes) == 5:
    #    arr=np.zeros(3)
    #    arr[principalaxes[0]] = eigenvalues[0]
    #    arr[principalaxes[1]] = eigenvalues[1]
    #    arr[principalaxes[2]] = eigenvalues[2]
    #    return arr
    #else:
    #    if principalaxes[0] == 2:
    #        return np.array([eigenvalues[1],eigenvalues[2],eigenvalues[0]])
    #    elif principalaxes[1] == 2:
    #        return np.array([eigenvalues[0],eigenvalues[2],eigenvalues[1]])
    #    elif principalaxes[2] == 2:
    #        return np.array([eigenvalues[0],eigenvalues[1],eigenvalues[2]])
    #    else:
    #        sys.exit("error, sort eigenstate is unavailable, no state having z principal axis")

def getCCEinDiagData(data,rxyzs,addspins,jt,Aio,Qio,Dio,Adata,Qdata,Ddata):
    # check addspin
    addspins=getSpinNames(addspins)

    # get A,Q data
    for i,atm in enumerate(rxyzs):
        for j,addspin in enumerate(addspins):
            val = re.findall('[a-zA-Z]+',addspin) # spin name
            if atm[0] == val[0]:
                data['rxyz'] = appendrxyz(data['rxyz'],jt,addspin,atm)
                if Aio:
                    A = setAtensor(Adata,addspin,i)
                    w,v = np.linalg.eig(A)
                    sort_w = sortEigXYZ_3x3(w,v)
                    data['Adiag'] = appendtensordiag(data['Adiag'],jt,addspin,sort_w)
                if Qio:
                    Q  = setQtensor(Qdata,i) 
                    w,v = np.linalg.eig(Q)
                    sort_w = sortEigXYZ_3x3(w,v)
                    data['Qdiag'] = appendtensordiag(data['Qdiag'],jt,addspin,sort_w)
                    data['eq'] = appendeq(data["eq"],addspin,Qdata['ceQ'][addspin])

    # get D data
    if Dio:
        D  = setDtensor(Ddata) 
        w,v = np.linalg.eig(D)
        sort_w = sortEigXYZ_3x3(w,v)
        data['Ddiag'] = appendtensordiag(data['Ddiag'],jt,'e',sort_w)

    return data

def readPreviousFile(data,fin):

    # read the cce.in if existing
    cceindata = readfile(fin,passio=True)

    for i,dat in enumerate(cceindata):
        tag=re.findall('[a-zA-Z]+',dat[0])[0]
        jt=int(re.findall('[0-9]+',dat[0])[0])
        vals=list(map(float,dat[3:]))
        if 'rxyz' == tag:
            data["rxyz"].append([jt,dat[2]]+vals)
        if 'A' == tag:
            data["A"].append([jt,dat[2]]+vals)
        if 'Q' == tag:
            data["Q"].append([jt,dat[2]]+vals)
        if 'ZFS' == tag:
            data["D"].append([jt,dat[2]]+vals)
        if 'Adiag' == tag:
            data["Adiag"].append([jt,dat[2]]+vals)
        if 'Qdiag' == tag:
            data["Qdiag"].append([jt,dat[2]]+vals)
        if 'Ddiag' == tag:
            data["Ddiag"].append([jt,dat[2]]+vals)

    return data

def writeNewFile(fout,data_ccein,data,megs,diag=False):

    # keywords that you want to add
    keys = ['rxyz','A','Q','D']
    if diag:
        keys = ['rxyz','Adiag','Qdiag','Ddiag']


    # write the data in cce.in
    f = open(fout,'w')
   
    for k in keys:
        meg=megs[k]
        for i,dat_ccein in enumerate(data_ccein[k]):
            f.write(meg.format(*dat_ccein)) 
        for j,dat in enumerate(data[k]):
            f.write(meg.format(*dat)) 

    f.close()

    
