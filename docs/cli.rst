Command Line Interface
======================

.. argparse::
   :module: crema.params
   :func: _configure_parser
   :prog: crema

Output
---------
Crema will produce five output files:

#. A text file named "crema.psms.txt" containing an additional column (crema q-value) of PSM-level confidence estimate results appended to the spectrum, peptide, and score columns from the input file.

#. A text file named "crema.peptides.txt" containing an additional column (crema q-value) of peptide-level confidence estimate results appended to the spectrum, peptide, and score columns from the input file.

#. A text file named "crema.proteins.txt" containing an additional column (crema q-value) of protein-level confidence estimate results appended to the protein ID and score columns from the input file.

#. A text file named "crema.protein_groups.txt" containing an additional column (crema q-value) of protein-level confidence estimate results appended to the protein ID and score columns from the input file.

#. A text file named "crema.log.txt" containing the logging information from using crema. Note that this file is only created when running crema through via the command line.

Note that the program writes to the current working directory by default.

* The name of the output directory can be specified using the ---output_dir argument.
* A prefix can be added to the output file names be using the ---file_root argument.
