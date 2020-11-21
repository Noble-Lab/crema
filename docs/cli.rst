Command Line Interface
======================

.. argparse::
   :module: crema.params
   :func: _configure_parser
   :prog: crema

Output
---------
There will be two output files:

#. A csv file named "crema.psm_results.txt" containing two columns (False Discovery Rate and Q-Value appended to the columns specified from the input file.
#. A csv file named "crema.logfile.log" containing the logging information from using crema at the CLI.

Note that the program writes to the current working directory by default:

* The name of the output directory can be specified using the ---output_dir argument.
* A prefix can be added to the output file name be using the ---file_root argument.