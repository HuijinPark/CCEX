# Getting Started
<a id="GettingStarted"></a>
  
This is the getting started guide.

<a id="HowToRun"></a>
## How to run

To run the code using <tt> mpirun </tt>, use the following command :  

``` 
$ mpirun -n process_number code_path -f input_file 
```  

- process_number : The number of processes to run.  
- code_path : The path to your executable code.  
- input_file : The input file required for the code to run. <br>
  (input file should be <tt>.json</tt> file)  

**Example** 

``` 
$ mpirun -n 4 /usr/local/bin/mycode -f data/input.json 
```

**Possble** **command-line** **options**
- -h, --help : Print the options and exit  
- -f : [Input file](InputFile.md) option (Required) (See [Input file](InputFile.md))  
- -m : Method (See [method](InputFile_General.md#tag_method))  
- -q : Quantity (See [quantity](InputFile_General.md#tag_quantity))   
- -I : File path to set the position of bath spins (See [bathfiles](InputFile_General.md#tag_bathfile))  
- -s : File path to set the spin state of bath spins (See [bath_state_file](InputFile_General.md#tag_bath_state_file))   
- -a : File path to set principal axies of the defects in the bath (See [Defect.axis_file](InputFile_Defect.md#tag_df_axis_file))  
- -S : File path to set the spin state of on-site spins of the defects in the bath (See [Defect.state_file](InputFile_Defect.md#tag_df_state_file))  
- -N : The number of state configuration to do time ensemble average (See [nstate](InputFile_General.md#tag_nstate))  
- -B : Magnetic field in z direction (Unit : G) (See [bfield](InputFile_General.md#tag_bfield))  
- -o : Result file path (See [outfile](InputFile_Output.md#tag_outfile))  
- -v : Verbose output (longer log)  
- -P : Write the cluster plan to `<outfile>.clusters.json` and stop before the dynamics (cCCE, gCCE and cCCE with `cluster_mode` pcce and explicit connections)  

Three run modes take the input file and replace the ordinary calculation:
- --auto_conv -f FILE : The automatic convergence; the same as `-f FILE` with an `auto_convergence` block (See [Automatic convergence tags](InputFile_AutoConvergence.md)). It takes no other option.  
- --auto_status -f FILE : Where that automatic convergence is (running, interrupted, stopped or how it ended); nothing is computed.  
- --adaptive ... : The adaptive refinement; every word after `--adaptive` goes to it (e.g. `--adaptive --input in.json --output DIR ...`). An `adaptive` block in the input asks for the same; giving both is refused (See [the `adaptive` block](../methods/adaptive-input-block.md)).  

**About** **options**

The code read options twice :
- Read options from input file called as In-file options
- Read optoins from command line called as Command-line options <br>
> **Note:** Command-line options are "priority"

<hr>
<a id="Workflow"></a>
## Workflow

The CCEX software follow the below work process. 

When you give the arguments in the command line, the code frist read the arguments. <br>
Note that if you don't give the input file accesible with -f option, then the code doesn't work. <br>

Code read the all tags in input file and update the ...  

<img src="CCEX_workflow.png" alt="" width="600">
