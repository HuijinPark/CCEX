#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from pyform import *
from textcolor import C
import os
import matplotlib.pyplot as pl
import numpy as np
from rwfiles import readfile, writefile
'''_________________________________________________________________'''

# Code Information
code_info = \
{
    'summary' : "Data Convergence test \n\n\
    This code include processes which are following:\n\
    1) Find Converged Point according to tolerance\n\
    2) Get coherence graph\n\
    3) Get difference graph\n\
    * Convergence test based on the following equation :\n\
      (1) reference data : coherence at the last condition (fully converged data)\n\
      (2) certain data : coherence at a certain condition (usually not converged data)\n\
      (3) area difference = sum(certain data - reference data)\n\
      (3) area of reference data = sum(reference data)\n\
      (4) difference(%) = (area difference) / (area of reference data) * 100\n\
      (5) Converged point : if difference < tolerance , then the variable would be converged point",
	'mkdate'  : f"{time.ctime(os.path.getctime(__file__))}", 
    'Author'  : "Huijin Park",
	'update1' : f"(2023.03.21) Key : {C.Y}--add{C.Z}\n\
    \t\t    If you want to average out the input data and convergence test for them\n\
    \t\t    Use this option", 
}
cal_form = cal_form(code_info['summary'],code_info['mkdate'])
cal_form._start(further_code_info=code_info)

'''_________________________________________________________________'''

# Input parameter Information
argv_info = \
{
'[1]' : '(s)fins          = Import file list',
'[2]' : '(s)variables     = variable range',
'[3]' : '(s)varialbe_name = variable name (e.g. rBath)',
'[4]' : '(s)unit          = unit',
'[5]' : '(f)tolerance     = Converged critria',
'[6]' : '(s)savefig       = outfile',
'[7]' : '(s,key)add       = Converged test with averaging out [ Opt. True / False{default} ]',
'[8]' : '(s,key)findconv  = Find converged data [ Opt. True{default} / False ]',
}
argv_add_meg = "\t\tInput file format\n\t\t : [t_free] [coherence_function] "

var = cal_form.is_argv_error(argv_info,argv_add_meg)

'''__________________________________________________________________'''


'''                               Main                               '''

'''__________________________________________________________________'''

fig, ax = pl.subplots(2,1,figsize=(6,8))

files = ((var.fins).strip()).split()
variables = list(map(float,((var.variables).strip()).split()))

if len(files) != len(variables):
    cal_form._end(meg_error='Error, lengthes of fins and variables are different')

# reference data
dat_ref = readfile(files[-1],ascomplex=True)

# Find difference
ISDONE=False
diffs=[]
THEVALUE=None
dat = None

# data would be added(prepare)
if var.add == 'True':
    dat=np.zeros(len(dat_ref[:,1]))

# find difference algorithm
for i,f in enumerate(files):
    
    # get data
    if var.add == 'True': 
        dat+=readfile(f,ascomplex=True)
    else:
        dat=readfile(f,ascomplex=True)
    variable=variables[i]
    ax[0].plot(dat[:,0],dat[:,1],label=f"{variable}{var.unit}",lw=2,ls=":")

    # data difference
    dat_masked = np.ma.masked_array(dat)
    dat_ref_masked = np.ma.masked_array(dat_ref)
    filled_array = dat_masked.filled(0) #masked array fill as 0
    filled_array_ref = dat_ref_masked.filled(0) #masked array fill as 0

    diff = np.sum(filled_array[:,1] - filled_array_ref[:,1])
    ref = np.sum(filled_array_ref[:,1])
    diffs.append(diff/ref*100)
    
    # Criteria
    print(variable, diffs[-1], abs(diffs[-1]))
    if abs(diffs[-1]) < var.tolerance and not ISDONE:
        THEVALUE=variable
        print(f"\t ! Converged_Point: {variable} ({var.unit})")
        if var.findconv != "False":
            ax[0].plot(dat[:,0],dat[:,1],label=f"{variable}{var.unit}",lw=2,ls="-",color="darkred")
        #ax[1].scatter(variable,diff,color="darkred",s=5)
        #ax[1].hlines(y=diff,xmin=variables[0],xmax=variable,ls=':',color='k',alpha=0.8)
        #ax[1].vlines(x=variable,ymin=0.0,ymax=diff,ls=':',color='k',alpha=0.8)
        ISDONE=True

ax[0].tick_params(axis='x',which='both',direction='in',pad=6,labelsize=20, top=True)
ax[0].tick_params(axis='y',which='both',direction='in',pad=6,labelsize=20, top=True)
ax[0].set_ylabel(r'$L (T)$',fontsize='23')
ax[0].set_xlabel(r'$T (ms)$',fontsize='23')
#ax[0].legend(ncol=3, loc='upper right',fancybox=False,prop={'size':20}, frameon=False)
ax[0].legend(ncol=1,bbox_to_anchor=(1.1,0.5),fancybox=False,prop={'size':20}, frameon=False)
#ax[0].set_xlim(0,0.3)

ax[1].vlines(x=THEVALUE,ymin=min(diffs),ymax=max(diffs),ls=':',color='k',alpha=0.8)
ax[1].plot(variables,diffs,lw=3,ls="-",markersize=5,marker='o',color="blue")
ax[1].tick_params(axis='x',which='both',direction='in',pad=6,labelsize=20, top=True)
ax[1].tick_params(axis='y',which='both',direction='in',pad=6,labelsize=20, top=True)
ax[1].set_ylabel(f'Difference (%)',fontsize='23')
ax[1].set_xlabel(f'{var.varialbe_name} ' + fr' ({var.unit})',fontsize='23')

pl.subplots_adjust(left=0.05, bottom=0.2, right=1, top=1.37, wspace=0.25, hspace=0.2)    
pl.savefig(var.savefig, dpi=300, facecolor='w', edgecolor='w',orientation='portrait', format=None,transparent=False, bbox_inches='tight', pad_inches=0.1)

cal_form._end()
