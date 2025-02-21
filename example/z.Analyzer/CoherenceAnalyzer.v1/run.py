import sys
import numpy as np
import copy
import logging
import json
from CoherencePlot_v2 import *
from DataSpace import *
import argparse

import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered in power")
warnings.filterwarnings("ignore", category=np.ComplexWarning, message="Casting complex values to real discards the imaginary part")
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

timescale={'ns':1e-9, 'us':1e-6, 'ms':1e-3,'s':1}

def get_parser():
    ####################################################################################################
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                    description="\
    ./DataSpace.py  -d  /home/huijin/cal/pycce/result/STRAIN/BFIELDG/rawdata/ \
                    -f  gCCE2_hBNVB_HE_BFIELD_conv10 \
                    -vn STRAIN BFIELD  \
                    -v  g95 g96 g97 flat t110  // Variable values for STRAIN \
                    -v  0 2 4 6 8  // Variable values for BFIELD")

    parser.add_argument('-d','--input_dir',dest='inputfile_dir',type=str,help="The directory of inputfile")
    parser.add_argument('-fi','--input',dest='inputfile',type=str,help="The inputfile name")
    parser.add_argument('-vn','--varnames',dest='variable_names',type=str,nargs='*',help="Variable names")
    parser.add_argument('-v','--vars',dest='variable_values',action='append',type=str,nargs='*',help="Variable values for each names")
    parser.add_argument('-ea','--ens_avg',dest='ensemble_average',action='store_true' ,help="Do ensemble average for all variables")
    parser.add_argument('-ntn','--nantozero',dest='nantozero',action='store_true' ,help="Do ensemble average for all variables")
    parser.add_argument('-uc','--unit_conv',dest='unit_conversion',default="ms", choices=['ms','us','ns','s'] ,help="Time unit conversion")

    #fitting parameter informations
    parser.add_argument('-fit_w','--fit_write',dest='write_fit',action='store_true', help="Write the fitting parameters")
    parser.add_argument('-pl_T','--pl_T2',dest='do_plot_T2',action='store_true', help="Plot for T2s")
    parser.add_argument('-pl_p','--pl_p',dest='do_plot_p',action='store_true', help="Plot for Powers")
    parser.add_argument('-fit_space','--fit_space',dest='fit_space',default=1,type=int, help="Plot for T2s")
    parser.add_argument('-fit_initp','--fit_initp',dest='fit_initp',default=None, type=float, help="Plot for Powers")

    # error correction
    parser.add_argument('-err_d','--err_corr_dir',dest='err_corr_dir',type=str,help="The directory of reference file")
    parser.add_argument('-err_fi','--err_corr_fi',dest='err_corr_fi',type=str,help="File name of the reference to check numerical error (noDiv)")
    parser.add_argument('-err','--err_corr',dest='err_corr', action='store_true', help="Error correation")
    parser.add_argument('-err_tol','--err_tol',dest='tolerance', default=2.0, type=float, help="Error correation tolerance")

    # plot
    parser.add_argument('-xlog','--xlog',dest='xlog',action='store_true', help="x axis log scale")
    parser.add_argument('-xlim','--xlim',dest='xlim',nargs='*', type=float, help="xlim = [xmin, xmax]")
    parser.add_argument('-ylim','--ylim',dest='ylim',nargs='*', type=float, help="ylim = [xmin, xmax]")
    parser.add_argument('-pl','--plot',dest='do_plot',action='store_true', help="Plot for all data")
    parser.add_argument('-pl_m','--plot_m',dest='do_plot_multi',action='store_true', help="Plot for all data in a figure")
    parser.add_argument('-fit','--fitting',dest='do_fitting',action='store_true', help="Show fitting on plot")
    parser.add_argument('-fo','--fout',dest='outfile',nargs='*',type=str, help="output file")

    # intensity plot
    parser.add_argument('-pl_i','--plot_intensity',dest='do_plot_intensity',action='store_true', help="Coherence Intensity Plot for all data ( x axis : time, y axis : variable_values )")
    parser.add_argument('-interpol_x','--interpolate_x',dest='interpolate_x',default=1, type=int, help="Interpolation of the plot (intensity plot only, usually time) ")
    parser.add_argument('-interpol_y','--interpolate_y',dest='interpolate_y',default=1, type=int, help="Interpolation of the plot (intensity plot only, usually b-field)")
   
    # verbosity
    parser.add_argument('-vb','--verbosity',dest='verbosity',action='store_true', help="Verbosity")

    return parser.parse_args()
    ####################################################################################################

def pretty_print_args(args):
    args_dict = vars(args)
    for key, value in args_dict.items():
        print(f"\t{key:20s}: {value}")
    print("")

def check_argument_error(args):
    try:
        if (len(args.variable_names) == len(args.variable_values)):
            pass;
        else:
            print("variable_names : ",args.variable_names,end='\n\n') 
            print("variable_values : ",args.variable_values,end='\n\n') 
            sys.exit("variable_values length is different from variable_names")
    
    except TypeError:
        print("variable_names : ",args.variable_names,end='\n\n') 
        print("variable_values : ",args.variable_values,end='\n\n') 
        sys.exit("TypeError")

def get_variables(args):
    variables = []
    for i,varname in enumerate(args.variable_names):
        variables.append([args.variable_names[i],args.variable_values[i]])
    
    return variables    


if __name__ == "__main__":

    ####################
    # Read arguments 
    ####################
    args = get_parser()
    print("\t","_"*50)
    print("\n\t< Input options >\n")
    pretty_print_args(args)
    check_argument_error(args)
    variables = get_variables(args)
    print("\t","_"*50)

    ####################
    # Read files & Data
    ####################
    print("\n\t< Read file >\n")
    Dsp = DataSpace() # Data Space
    Dsp.verbosity = args.verbosity
    Dsp.get_filenames_format(args.inputfile_dir,args.inputfile) # Get file format
    Dsp.get_variable_all(variables) # Get variables for each variable names
    Dsp.update_filenames() # Update the file names with the variables

    count = Dsp.read_file_coherence_complex_v2(timeunit_converter=timescale['ms']/timescale[args.unit_conversion]) # Read the files
    print(f"\n\tRead file number = {count} ...")
    print("\t","_"*50)

    ####################
    # Nan to zero
    ####################
    if (args.nantozero):
        print("\n\t< Nan to Zero >\n")
        Dsp.refine_nan2num_coherence_all()
        print("\tDone to change the nan number to zero")
        print("\t","_"*50)

    ####################
    # Error correction 
    ####################
    Dsp_ref = DataSpace() # Reference data (no divided data)
    if (args.err_corr):
        print("\n\t< Error correction >\n")
        print("\tCorrect the numerical error ... \n")
        Dsp_ref.get_filenames_format(args.err_corr_dir,args.err_corr_fi)
        Dsp_ref.get_variable_all(variables)
        Dsp_ref.update_filenames()
        Dsp_ref.read_file_coherence_complex_v2(timeunit_converter=timescale['ms']/timescale[args.unit_conversion])

        # Do error correction
        for key,data in Dsp.data_dict.items():
            if not np.all(data[:,0] == Dsp_ref.data_dict[key][:,0]): sys.exit("Error(Error correction), the time scale is not the same")
            ts = data[:,0]
            L_wD = data[:,1]
            L_nD = Dsp_ref.data_dict[key][:,1]
            L_EC = Dsp._error_correction_coherence(L_wD,L_nD,args.tolerance) # do error correction
            Dsp.data_dict[key][:,1] = L_EC # change initial data into error corrected data

            # Write
            write_vars = Dsp._get_variables_from_key(key)
            fname = args.outfile[0]
            for i,write_varname in enumerate(Dsp.variables_keys):
                fname = fname.replace(write_varname,write_vars[i])

            print(f"\t\tWrite file : {fname}")
            writefile(fname,['f','c'],[ts.real,L_EC])
        print("\n\tDone to correct error")
        print("\t","_"*50)

        
    ####################
    # Ensemble average
    ####################
    ens_count = None 
    if (args.ensemble_average):
        print("\n\t< Ensemble average >\n")
        # Write
        ts, L_ens, ens_count = Dsp.average_coherence_all()
        fname = f"{args.outfile[0]}_conv{ens_count}"
        writefile(fname,['f','c'],[ts.real,L_ens])
        print(f"\n\t\tTotal averaged number : {ens_count} #")
        print(f"\t\tWrite ensemble averaged result : {fname}") 
        
        ts_2d = ts[:,np.newaxis]
        L_ens_2d = L_ens[:,np.newaxis]

        Dsp.data_dict["avg"] = np.concatenate((ts_2d,L_ens_2d),axis=1) 
        print("\n\tDone ensemble average")
        print("\t","_"*50)
   
    ####################
    # Single plot
    ####################
    if (args.do_plot):
        print("\n\t< Single Plot >\n")

        print(f"\t\ttime unit : {args.unit_conversion}")
        # show with fitting
        if args.do_fitting:
            print(f"\t\tShow fitting (space : {args.fit_space} , initial_power : {args.fit_initp})")

        # do single plot
        if (len(args.variable_names) == 1) and (not args.ensemble_average):
            print("\n\t\tPlot each data in each figure")
            for i,varname in enumerate(args.variable_names):
                for j,varval in enumerate(args.variable_values[i]):
                    Fig = Dsp.plot_coherence_one(varname,varval,{} ,variable_name=varval,xunit=args.unit_conversion,colors='k',showWfitting=args.do_fitting,space=args.fit_space,init_fitPower=args.fit_initp)

                    if (args.xlim != None):
                        Fig.Fig.ax.set_xlim(args.xlim[0],args.xlim[1])
                    if (args.ylim != None):
                        Fig.Fig.ax.set_ylim(args.ylim[0],args.ylim[1])

                    if (args.xlog):
                        Fig.Fig.ax.set_xscale('log')

                    fname=f"{args.outfile[0]}".replace(varname,varval)
                    Fig.Fig.save(fname)
                    if args.verbosity:
                        print(f"\n\t\tSave the plot in {fname}")

        # do ensemble averaged single plot 
        elif  args.ensemble_average:
            print("\t\tPlot ensemble averaged data")
            varname="None"
            varval="avg"
            Fig = Dsp.plot_coherence_one(varname,varval,{} ,variable_name=varval,xunit=args.unit_conversion,showWfitting=args.do_fitting,space=args.fit_space,init_fitPower=args.fit_initp)

            if (args.xlim != None):
                Fig.Fig.ax.set_xlim(args.xlim[0],args.xlim[1])
            if (args.ylim != None):
                Fig.Fig.ax.set_ylim(args.ylim[0],args.ylim[1])
            if (args.xlog):
                Fig.Fig.ax.set_xscale('log')

            fname=f"{args.outfile[0]}_conv{count}"
            Fig.Fig.save(fname)
            print(f"\n\t\tSave the plot in {fname}")

        else: 
            sys.exit("Error args.variable_names is larger than 1")
        print("\n\tDone single plot")
        print("\t","_"*50)

    ####################
    # Multi plot
    ####################
    if (args.do_plot_multi):
        print("\n\t< Multi Plot >\n")
        print(f"\t\ttime unit : {args.unit_conversion}")

        if len(args.variable_names) == 1:
            colors=get_cmap_colors(len(args.variable_values[0]))
            varname=args.variable_names[0]
            Dsp.plot_coherence_all(varname,{},variable_name=args.variable_values[0],xunit=args.unit_conversion,colors=colors)
            Fig = Dsp.fig_coherence
            if (args.xlim != None):
                Fig.Fig.ax.set_xlim(args.xlim[0],args.xlim[1])
            if (args.ylim != None):
                Fig.Fig.ax.set_ylim(args.ylim[0],args.ylim[1])

            if (args.xlog):
                Fig.Fig.ax.set_xscale('log')
            
            Fig.Fig.legend(0.9,0.9,frameon=False)
            Fig.Fig.save(f"{args.outfile[0]}")
        else: sys.exit("Error args.variable_names is larger than 1")
        print("\n\tDone Multi plot")
        print("\t","_"*50)

    ####################
    # Intensity plot 
    ####################
    if (args.do_plot_intensity):
        print("\n\t< Intensity Plot >\n")
        print(f"\t\ttime unit : {args.unit_conversion}")
        print(f"\t\tinterpolation : [{args.interpolate_x} , {args.interpolate_y}]")
        if len(args.variable_names) == 1:
            for i, varname in enumerate(args.variable_names):
                varlist = Dsp.variables_dict[varname]
                tlist = Dsp.data_dict[f'{varlist[0]}'][:,0]
                Condition = {}
                interpolate = [args.interpolate_x,args.interpolate_y]
                Iplot = Dsp.plot_coherence_intensity(varname,varlist,Condition,cmap="RdYlBu_r",xunit=args.unit_conversion, interpolate=(len(varlist)*interpolate[0],len(tlist)*interpolate[1]))
                if (args.xlim != None):
                    Iplot.ax.set_xlim(args.xlim[0],args.xlim[1])
                if (args.ylim != None):
                    Iplot.ax.set_ylim(args.ylim[0],args.ylim[1])


                Iplot.save(f"{args.outfile[0]}")
        else: sys.exit("Error args.variable_names is larger than 1")
        print("\n\tDone intensity plot")
        print("\t","_"*50)

    ####################
    # Write fitting parameters 
    ####################
    if (args.write_fit):
        print("\n\t< Wrtie the fitting result >\n")
        print(f"\t\ttime unit : {args.unit_conversion}")
        print(f"\t\tFitting with space : {args.fit_space} , initial_power : {args.fit_initp}")
        print("\t\tFormat : variable in order, T2, power\n")
        if len(args.variable_names) == 1:
            Dsp.do_fitting_coherence_complex(space=args.fit_space, init_fitPower=args.fit_initp)

            t2s = list(Dsp.t2_dict.values())
            ps  = list(Dsp.power_dict.values())
            var = np.asarray(args.variable_values[0]).astype(float)

            # Check none
            t2s_new = []
            ps_new = []
            var_new = []
            for i, t2 in enumerate(t2s):
                if t2!=None:
                    t2s_new.append(t2s[i])
                    ps_new.append(ps[i])
                    var_new.append(var[i])
                else:
                    print(f"\n\t\t! When {args.variable_names[0]} is {var[i]}, the fitting cannot be done")

            params = list(zip(t2s_new, ps_new, var_new))
            
            params_sorted = sorted(params, key=lambda x: x[2])

            t2s_sorted = [x[0] for x in params_sorted]
            ps_sorted = [x[1] for x in params_sorted]
            var_sorted = [x[2] for x in params_sorted]

            writefile(f"{args.outfile[0]}_params",['f','f','f'],[var_sorted,t2s_sorted,ps_sorted])

            # plot
            fig = EasyPlot()
            if args.do_plot_T2:
                fig.marker='o'
                fig.ls='--'
                fig.plot(var_sorted,t2s_sorted)
                if (args.xlim != None):
                    fig.ax.set_xlim(args.xlim[0],args.xlim[1])
                if (args.ylim != None):
                    fig.ax.set_ylim(args.ylim[0],args.ylim[1])
                fig.save(f"{args.outfile[0]}_T2")

            fig = EasyPlot()
            if args.do_plot_p:
                fig.marker='o'
                fig.ls='--'
                fig.plot(var_sorted,ps_sorted)
                if (args.xlim != None):
                    fig.ax.set_xlim(args.xlim[0],args.xlim[1])
                if (args.ylim != None):
                    fig.ax.set_ylim(args.ylim[0],args.ylim[1])
                fig.save(f"{args.outfile[0]}_Power")

        else: sys.exit("Error args.variable_names is larger than 1")
        print("\n\tDone to fit")
        print("\t","_"*50)

        
            
