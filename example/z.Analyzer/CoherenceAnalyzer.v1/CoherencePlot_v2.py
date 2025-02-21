# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 15:26:52 2023

@author: 17732
"""


import matplotlib as mpl
import matplotlib.pyplot as pl
import numpy as np
from scipy.interpolate import interp2d, interp1d
import copy

def show_possible_cmap():
    print('Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r')
    
def get_cmap(cmap='jet'):
    return pl.cm.get_cmap(cmap)

def get_cmap_color(totaln,i,cmap='jet'):
    return pl.cm.get_cmap(cmap)(np.linspace(0,1,totaln))[i]

def get_cmap_colors(totaln,cmap='jet'):
    return pl.cm.get_cmap(cmap)(np.linspace(0,1,totaln))
            
class EasyPlot:
    
    def __init__(self,figsize = None,nrow=1,ncol=1):
        
        self.figsize = figsize
        
        if self.figsize != None:
            self.fig = pl.figure(figsize=self.figsize) #
        else:
            self.fig = pl.figure(figsize=(5,4)) #
        
        self.nrow = nrow
        self.ncol = ncol
        self.ax_i = 1
        self.ax  = (self.fig).add_subplot(self.nrow,self.ncol,self.ax_i)
        self.axes = [self.ax]
    
        self.marker=""
        self.ls="-"
        self.markersize=2
        self.lw=1.5
        self.alpha=0.9
        self.label=""
        self.pad=6
        self.labelsize=15
        self.xlabel=""
        self.ylabel=""
        self.fontsize=20
        self.xlim=[]
        self.ylim=[]
        self.bbox_to_anchor=() #(0.0,1.0)
        self.color=None
        # self.cmap=pl.cm.get_cmap("Set2")
        self.title=""

            
    def plot(self,x,y,background=False,legend=False,ax_i=1):
        
        if len(self.axes)!=1 or ax_i!=1:
            try:
                self.axes[ax_i-1]
                self.ax_i = ax_i
                self.ax = self.axes[ax_i-1]
            except:
                self.ax_i = ax_i
                self.ax = (self.fig).add_subplot(self.nrow,self.ncol,self.ax_i)
                self.axes.append(self.ax)
        
        line = ""
        try:
            if self.color == None:
                line = (self.ax).plot(x,y,marker=self.marker,markersize=self.markersize,ls=self.ls,lw=self.lw,alpha=self.alpha,label=self.label)
            else:
                line = (self.ax).plot(x,y,marker=self.marker,markersize=self.markersize,ls=self.ls,lw=self.lw,alpha=self.alpha,label=self.label,color=self.color)
        except:
            line = (self.ax).plot(x,y,marker=self.marker,markersize=self.markersize,ls=self.ls,lw=self.lw,alpha=self.alpha,label=self.label,color=self.color)
        
        (self.ax).tick_params(axis='both',which='both',direction='in',pad=self.pad,labelsize=self.labelsize, top=True)
        (self.ax).set_xlabel(self.xlabel,fontsize=self.fontsize,fontname="Arial")
        (self.ax).set_ylabel(self.ylabel,fontsize=self.fontsize,fontname="Arial")
        (self.ax).set_title(self.title,fontsize=self.fontsize,fontname="Arial")
        
        if background:
            (self.ax).grid(b=True, which='major', color='#666666', linestyle='--',alpha=0.5)
            (self.ax).minorticks_on()
            (self.ax).grid(b=True, which='minor', color='#999999', linestyle='--', alpha=0.2)
        
        if self.xlim != []:
            (self.ax).set_xlim(self.xlim[0],self.xlim[1])
            
        if self.ylim != []:
            (self.ax).set_ylim(self.ylim[0],self.ylim[1])
            
        # self.ax.tick_params(axis='x',which='major',direction='in',pad=self.pad,labelsize=self.labelsize, top=True)
        
        return line
    
    def legend(self,x,y,fontsize=12,frameon=True):
        self.ax.legend(fontsize=fontsize,bbox_to_anchor=(x,y),frameon=frameon)
        
    def do_rcParams(self,key,value):
        mpl.rcParams[key]=value
        
    def show_possible_rcParams(self):
        # print(mpl.rcParams)
        for k,v in mpl.rcParams.items():
            print(f"{k:30s} : {v}")
        
        print("Usually used")
        print("   legend.framealpha ")
        print("   legend.facecolor ")
        print("   legend.edgecolor ")
        print("   legend.fontsize ")
        print("   legend.loc ")
        print("   axes.edgecolor ")
        print("   xtick.bottom ")
        print("   xtick.left ")
        print("   xtick.labelsize ")
        
   
    def update_fig(self):
        self.fig.canvas.draw()
        
    def save(self,fout):
        #print(f"Plot save : {fout}")
        #self.fig.tight_layout()
        self.fig.savefig(fout,dpi=500, bbox_inches='tight')
 
    def scatter(self,x,y,background=False):
        
        if self.color == None:
            (self.ax).scatter(x,y,marker=self.marker,s=self.markersize,alpha=self.alpha,label=self.label)
        else:
            (self.ax).scatter(x,y,marker=self.marker,s=self.markersize,alpha=self.alpha,label=self.label,color=self.color)
   
        (self.ax).tick_params(axis='both',which='both',direction='in',pad=self.pad,labelsize=self.labelsize, top=True)
        (self.ax).set_xlabel(self.xlabel,fontsize=self.fontsize,fontname="Arial")
        (self.ax).set_ylabel(self.ylabel,fontsize=self.fontsize,fontname="Arial")
        (self.ax).set_title(self.title,fontsize=self.fontsize,fontname="Arial")
        
        if background:
            (self.ax).grid(b=True, which='major', color='#666666', linestyle='--',alpha=0.5)
            (self.ax).minorticks_on()
            (self.ax).grid(b=True, which='minor', color='#999999', linestyle='--', alpha=0.2)
        
        if self.xlim != []:
            (self.ax).set_xlim(self.xlim[0],self.xlim[1])
            
        if self.ylim != []:
            (self.ax).set_ylim(self.ylim[0],self.ylim[1])
        
    def show_possible_cmap(self):
        print('Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r')
        
    def get_cmap(self,cmap='jet'):
        return pl.cm.get_cmap(cmap)
    
    def get_cmap_color(self,totaln,i,cmap='jet'):
        return pl.cm.get_cmap(cmap)(np.linspace(0,1,totaln))[i]
        
class IntensityPlot:
    
    
    def __init__(self,figsize = None):
        
        self.figsize = figsize
        
        if self.figsize != None:
            self.fig = pl.figure(figsize=self.figsize) #
        else:
            self.fig = pl.figure(figsize=(7,5)) #
            
        self.ax  = (self.fig).add_subplot(1,1,1)
        
        self.marker=""
        self.ls="-"
        self.markersize=2
        self.lw=1.5
        self.alpha=0.9
        self.label=""
        self.pad=6
        self.labelsize=15
        self.fontsize=20
        self.xlim=[]
        self.ylim=[]
        self.clim=[0,1]
        self.bbox_to_anchor=() #(0.0,1.0)
        self.cmap='RdYlBu_r'
        self.title=""
        self.nfig=0
    
    def intensityplot(self,t,B,Lset,ylabel="B",xunit=r'$\mu$s',yunit="G"):
        
        X, Y = np.meshgrid(t, B)
        Z = np.zeros((len(B),len(t)))
        
        pcm  = self.ax.pcolormesh(t, B, Lset, cmap=self.cmap,vmin=self.clim[0], vmax=self.clim[1]);
        (self.ax).set_xlabel(rf"2$\tau$({xunit})",fontsize=self.fontsize,fontname="Arial")
        (self.ax).set_ylabel(rf"{ylabel}({yunit})",fontsize=self.fontsize,fontname="Arial")
        (self.ax).set_title(self.title,fontsize=self.fontsize,fontname="Arial")
        (self.ax).tick_params(axis='both',which='both',direction='in',pad=self.pad,labelsize=self.labelsize, top=True)
        cbar = self.fig.colorbar(pcm, ax=self.ax)
        cbar.set_label('Re[L]',fontsize=18)
                
        return pcm

    def interpolate_array(self,input_array, new_shape):
        
        if input_array.shape == new_shape:
            return input_array
        n, m = input_array.shape
        x = np.linspace(0, m - 1, m)
        y = np.linspace(0, n - 1, n)
        f = interp2d(x, y, input_array, kind='cubic')

        new_x = np.linspace(0, m - 1, new_shape[1])
        new_y = np.linspace(0, n - 1, new_shape[0])
        interpolated_array = f(new_x, new_y)

        return interpolated_array

    def interpolate_1d_array(self,input_array, new_length):
        
        if input_array.shape[0] == new_length:
            return input_array
        
        n = len(input_array)
        x = np.linspace(0, n - 1, n)
        f = interp1d(x, input_array, kind='cubic')

        new_x = np.linspace(0, n - 1, new_length)
        interpolated_array = f(new_x)

        return interpolated_array

    def save(self,fout):
        self.fig.tight_layout(pad=2.0)
        self.fig.savefig(fout,dpi=500)

##############################################################################

class EasyPlotCoherence:
    
    def __init__(self,unit=r"$\mu$s",figsize = None):
        
        self.Fig = EasyPlot(figsize=figsize)
        self.Fig.ls="-"
        self.Fig.lw='1.5'
        self.Fig.xlabel = fr"t({unit})"
        self.Fig.ylabel = "Coherence"
        self.Fig.alpha = 0.5
        self.unit = unit
   
    def do(self,t,L,variable=""):

        if variable!="":
            self.Fig.label=variable
        line = self.Fig.plot(t,L)

        self.fig = self.Fig.fig
        self.ax  = self.Fig.ax
        self.t = t
        self.L = L
        return line
    
    def scatter(self,t,L,variable=""):

        if variable!="":
            self.Fig.label=variable
        self.Fig.scatter(t,L)

        # self.fig = self.Fig.fig
        # self.ax  = self.Fig.ax
        # self.t = t
        # self.L = L
    
    def legend(self,x,y,frameon=True):
        # self.Fig.ax.legend(fontsize=12,bbox_to_anchor=(x,y))
        self.Fig.ax.legend(fontsize=12,loc=(x,y),frameon=frameon)
        
    def show_params(self,T2,n,showT2=True,showPower=True,adjust_x=0.0,adjust_y=0.0,adjust_space=0.0):

        t = self.t
        L = self.L
        
        if showT2:
            self.Fig.ax.text(t[-1]*0.05+adjust_x,L[0]*0.3+adjust_y,fr"$T_2$ = {T2:.4f}{self.unit}",fontsize=self.Fig.pad*2)
        if showPower:
            self.Fig.ax.text(t[-1]*0.05+adjust_x,L[0]*0.2+adjust_y+adjust_space,fr"$n$ = {n:.2f}",fontsize=self.Fig.pad*2)
            
    def get_cmap(self,cmap='jet'):
        return self.Fig.get_cmap(cmap)

    def show_possible_cmap(self):
        self.Fig.show_possible_cmap()


class ErrorBarPlot:
    
    def __init__(self,xlabel=r"$B_0(G)$",ylabel = r"$T_{2}(\mu s)$",figsize = None,nrow=1,ncol=1):   
        
        self.figsize = figsize
        
        if self.figsize != None:
            self.fig = pl.figure(figsize=self.figsize) #
        else:
            self.fig = pl.figure(figsize=(7,5)) #
        
        self.nrow = nrow
        self.ncol = ncol
        self.ax_i = 1
        self.ax  = (self.fig).add_subplot(self.nrow,self.ncol,self.ax_i,constrained_layout=True)
        self.axes = [self.ax]
        
        self.marker=""
        self.ls="-"
        self.markersize=2
        self.lw=1.5
        self.alpha=0.9
        self.label=""
        self.pad=6
        self.labelsize=15
        self.xlabel=""
        self.ylabel=""
        self.fontsize=20

        self.bbox_to_anchor=() 
        self.color=None
        self.title=""
        self.yerr=0.0
        self.xerr=0.0
        self.capsize=5
        self.ecolor=None

   
    def do(self,x,y,variable="",ax_i=1):
         
        if len(self.axes)!=1 or ax_i!=1:
            try:
                self.axes[ax_i-1]
                self.ax_i = ax_i
                self.ax = self.axes[ax_i-1]
                        
            except:
                self.ax_i = ax_i
                self.ax = (self.fig).add_subplot(self.nrow,self.ncol,self.ax_i)
                self.axes.append(self.ax)   

        try:
            if self.color ==None:
                self.ecolor='r'
        except:
            self.ecolor=copy.copy(self.color)     

        if variable!="":
            try:
                self.ax.errorbar(x,y,yerr=self.yerr,xerr=self.xerr,ls=self.ls,alpha=self.alpha,color=self.color,ecolor=self.ecolor,capsize=self.capsize,label=variable,marker=self.marker,markersize=self.markersize,lw=self.lw)
            except:
                self.ax.errorbar(x,y,yerr=self.yerr,xerr=self.xerr,ls=self.ls,alpha=self.alpha,ecolor=self.ecolor,capsize=self.capsize,label=variable,marker=self.marker,markersize=self.markersize,lw=self.lw)
        else:
            try:
                self.ax.errorbar(x,y,yerr=self.yerr,xerr=self.xerr,ls=self.ls,alpha=self.alpha,color=self.color,ecolor=self.ecolor,capsize=self.capsize,marker=self.marker,markersize=self.markersize,lw=self.lw)
            except:
                self.ax.errorbar(x,y,yerr=self.yerr,xerr=self.xerr,ls=self.ls,alpha=self.alpha,ecolor=self.ecolor,capsize=self.capsize,marker=self.marker,markersize=self.markersize,lw=self.lw)
                
        self.ax.set_xlabel(self.xlabel,fontsize=self.fontsize,fontname="Arial")
        self.ax.set_ylabel(self.ylabel,fontsize=self.fontsize,fontname="Arial")
        (self.ax).tick_params(axis='both',which='both',direction='in',pad=self.pad,labelsize=self.labelsize)
        
    
    def legend(self,x,y,frameon=True):
        # self.Fig.ax.legend(fontsize=12,bbox_to_anchor=(x,y))
        self.ax.legend(fontsize=12,loc=(x,y),frameon=frameon)

    
    def save(self,fout):
        # print(fout)
        self.fig.savefig(fout,dpi=500, bbox_inches='tight')
 

class EasyPlotCoherence_T2:
    
    def __init__(self,xlabel,xunit="",yunit=r"$\mu$s",ylabel=r"$T_2$",figsize = None):
        
        self.Fig = EasyPlot(figsize=figsize)
        self.Fig.xlabel = f"{xlabel}({xunit})"
        self.Fig.ylabel = fr"{ylabel}({yunit})"
        self.Fig.marker="o"
        self.Fig.ls=":"
        self.Fig.markersize=4
    
    def do(self,xvars,t2s,variable=""):
        
        if variable!="":
            self.Fig.label=variable
        line = self.Fig.plot(xvars,t2s)
            
        self.fig = self.Fig.fig
        self.ax  = self.Fig.ax    
        return line


class EasyPlotCoherence_Power:
    
    def __init__(self,xlabel,xunit="",yunit=r"$\mu$s",figsize = None):
        
        self.Fig = EasyPlot(figsize=figsize)
        self.Fig.xlabel = f"{xlabel}({xunit})"
        self.Fig.ylabel = fr"$Power$"
        self.Fig.marker="o"
        self.Fig.ls=":"
        self.Fig.markersize=4
    
    def do(self,xvars,powers,variable=""):
        
        if variable!="":
            self.Fig.label=variable
        line = self.Fig.plot(xvars,powers)
            
        self.fig = self.Fig.fig
        self.ax  = self.Fig.ax  
        
        return line

class Histogram:

    def __init__(self,xlabel=r"$T_{2}(\mu s)$",ylabel = "Frequency(#)",figsize = None):   
        
        self.figsize = figsize
        
        if self.figsize != None:
            self.fig = pl.figure(figsize=self.figsize) #
        else:
            self.fig = pl.figure(figsize=(7,5)) #
            
        self.ax  = (self.fig).add_subplot(1,1,1)

        self.alpha=0.6
        self.label=""
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.labelsize=14
        self.fontsize=16
        self.pad=8
        self.xlim=[]
        self.ylim=[]
        self.bbox_to_anchor=() #(0.0,1.0)
        self.color='gray'
        self.title=""
        self.nfig=0
    
    def do(self,vals,bins,density=True):
        self.ax.hist(vals, bins=bins, density=density, alpha=self.alpha, edgecolor='k', color=self.color)
        self.ax.set_xlabel(self.xlabel,fontsize=self.fontsize,fontname="Arial")
        self.ax.set_ylabel(self.ylabel,fontsize=self.fontsize,fontname="Arial")
        (self.ax).tick_params(axis='both',which='both',direction='in',pad=self.pad,labelsize=self.labelsize)
