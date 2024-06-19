%CCEX-1.0.0 {#mainpage}
=================

Generated on 2024.06.01 

Introduction
============

This page provides a way how to simulate the spin qubit dynamics with the CCEX software.
The source code of this software can be found in [here](https://github.com/HuijinPark/CCEX).

About the simulation
============

In this code, we simulate the qubit dynamics in the presence of the bath spins.

You can simulate the spin qubit dynamcis with the following methods:
- Cluster correlation expension method (CCE method)
- CCE-NeNn method
- Generalized CCE method
- Partition CCE method (To be done)

You can obtain the following physical quantities of qubit for evolving time:
- Coherence function
- Density matrix (gCCE method only)
- Noise function (To be done for CCE, gCCE)

Installation
============

```
$ git clone https://github.com/HuijinPark/CCEX.git
```

## Requirement
[Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)
[Intel mpi](https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html#gs.b1bib5)

## Note
MakeFile Path should be changed.. (To be changed)

Usage Instructions
============

To run the code using `mpirun`, use the following command :

```
$ mpirun -n process_number code_path -f input_file

```

- process_number : The number of processes to run.
- code_path : The path to your executable code.
- \-f input_file : The input file required for the code to run. 
   * (input file should be `.json`<br> file)

# example 

```
$ mpirun -n 4 /usr/local/bin/mycode -f data/input.json
```

Update
============

`To be written.` 

To do list
===========
- Grouping of files in Model / Parser / Generator categories
- Index files based on IndexIntf
  - HTML navigation
  - HTML Help (chm)
  - Documentation Sets (XCode)
  - Qt Help (qhp)
  - Eclipse Help
- Search index
  - JavaScript based
  - Server based
  - External
- Citations
  - via bibtex
- Various processing steps for a comment block
  - comment conversion
  - comment scanner
  - markdown processor
  - doc tokenizer
  - doc parser
  - doc visitors
- Diagrams and Images
  - builtin
  - via Graphviz dot
  - via mscgen
  - PNG generation
- Output formats: OutputGenerator, OutputList, and DocVisitor
  - Html:  HtmlGenerator and HtmlDocVisitor
  - Latex: LatexGenerator and LatexDocVisitor
  - RTF:   RTFGenerator and RTFDocVisitor
  - Man:   ManGenerator and ManDocVisitor
  - XML:   generateXML() and XmlDocVisitor
  - print: debugging via PrintDocVisitor
  - text:  TextDocVisitor for tooltips
  - perlmod
- i18n via Translator and language.cpp
- Customizing the layout via LayoutDocManager
- Parsers
  - C Preprocessing
    - const expression evaluation
  - C link languages
  - Python
  - Fortran
  - VHDL
  - Tag files
- Marshaling to/from disk
- Portability functions
- Utility functions

