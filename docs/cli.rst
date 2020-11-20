Command Line Interface
======================

.. argparse::
   :module: crema.params
   :func: _configure_parser
   :prog: crema

Output
---------
A csv file named "crema.psm_results.txt" containing two columns (False Discovery Rate and Q-Value)
appended to the columns specified from the input file.

Note that the program writes to the current working directory by default.
The name of the output directory can be specified using the ---output_dir argument.
A prefix can be added to the output file name be using the ---file_root argument.