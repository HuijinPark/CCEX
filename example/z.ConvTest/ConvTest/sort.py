from rwfiles import readfile, writefile
import sys

fin=sys.argv[1]

dat=readfile(fin)

for i in range(len(dat)):
    dat[i][0] = int(dat[i][0])

outdat = sorted(dat,key=lambda l:l[0])

x=[]
y=[]
z=[]
w=[]
for i in range(len(outdat)):
    x.append(outdat[i][0])
    y.append(outdat[i][1])
    z.append(outdat[i][2])
    w.append(outdat[i][3])

writefile(fin,['d','s','s','s'],[x,y,z,w],)
