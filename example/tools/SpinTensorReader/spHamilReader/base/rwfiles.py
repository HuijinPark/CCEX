import spHamilReader.base.pyform as pyform
import numpy as np
'''_________________________________________________________________

    'summary' :Read any files"
   _________________________________________________________________'''

def readfile(fname,passio=None,prt=None,asfloat=None,asbath=None):
    data = []
    try: 
        open(fname, 'r')
    except FileNotFoundError:
        if passio !=None:
            pass;
        else: 
            pyform._error("\"%s\" file not found ! "%fname)
    else:
        with open(fname, 'r') as f:
            for i, line in enumerate(f):
                data.append(line.strip().split())

    if prt !=None:
        for i in range(len(data)):
            print(data[i])

    if asfloat!=None:
        return (np.asarray(data)).astype(float)

    if asbath!=None:
        dataArr = np.zeros((len(data),3))
        dataArr2 = []
        for i in range(len(data)):
            dataArr[i][0] = data[i][0]
            dataArr[i][1] = data[i][1]
            dataArr[i][2] = data[i][2]
            dataArr2.append(data[i][3])
        return (dataArr,dataArr2)

    return data



'''_________________________________________________________________

    'summary' :Write any files"
   _________________________________________________________________'''


def writefile(fout,dtype,data,comment=None):
    ''' form : 
                data  = [xArr,yArr,zArr, ... ] or list 
                dtype = [f,s,d, ... ]
                len(data) = len(dtype)
    '''
    f = open(fout, 'w')

    for col in range(len(data[0])):
        for row in range(len(data)):
            if dtype[row] == 'f':
                line = "{:>10.5f}\t".format(data[row][col])
            elif dtype[row] == 'd':
                line = "{:>10d}\t".format(data[row][col])
            elif dtype[row] == 's':
                line = "{:>10s}\t".format(data[row][col])
            f.write(line)
        line = "\n"
        f.write(line)
    
    if comment !=None:
        f.write(comment)

    f.close()

