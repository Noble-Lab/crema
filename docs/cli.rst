Command Line Interface
======================

.. argparse::
   :module: crema.params
   :func: _configure_parser
   :prog: crema

Output
---------
There will be three output files:

#. A txt file named "crema.psms.txt" containing an additional column (crema q-value) of psm level confidence estimate results appended to the spectrum, peptide, and score columns from the input file.

#. A txt file named "crema.peptides.txt" containing an additional column (crema q-value) of peptide level confidence estimate results appended to the spectrum, peptide, and score columns from the input file.

#. A txt file named "crema.log.txt" containing the logging information from using crema.

Note that the program writes to the current working directory by default:

* The name of the output directory can be specified using the ---output_dir argument.
* A prefix can be added to the output file names be using the ---file_root argument.