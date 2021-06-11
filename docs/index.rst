.. crema documentation master file, created by
   sphinx-quickstart on Fri Sep 11 15:54:39 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :hidden:
   :maxdepth: 2
   :titlesonly:
   :caption: crema

   self
   cli.rst
   api/index.rst
   Changelog <CHANGELOG.md>


Getting Started
---------------
**crema** produces confidence estimates for peptide detection in mass spectrometry proteomics experiments.
It takes as input files holding data regarding peptide-spectrum matches (PSMs), executes the
desired estimation method, and produces as ouptut confidence estimates of the PSMs.

Introduction
------------
One of the fundamental tasks in mass spectrometry proteomics is detecting peptides on the basis of the observed mass
spectra. Many tools exist to assign peptides to spectra, but unfortunately this matching is never 100% accurate,
meaning that there is uncertainty about whether a given PSM is correct or a false positive. We want to be able to
quantify this uncertainty so that we can be confident in our conclusions and ensure that expensive
downstream validation experiments use relevant and accurate data.

crema is a Python package that implements various methods to estimate false discovery rates (FDR)
in mass spectrometry proteomics experiments. crema focuses on
methods that rely on the concept of "target-decoy competition." The sole purposes of crema is to do decoy-based FDR
estimation, and to do it well. As a result, crema is lightweight and flexible. It has minimal dependencies and
supports a wide range of input and output formats. On top of that, it is extremely simple to use.

Ready to try crema for your analyses? See below for details on how to install and use crema.

Installation
------------
Before you can install and use crema, you'll need to have Python 3.6+
installed. If you think it may be installed, you can check with:

.. code-block:: bash

   $ python3 --version

If you need to install Python, we recommend using the `Anaconda Python
distribution <https://www.anaconda.com/products/individual>`_. This distribution
comes with most of the crema dependencies installed and provides the conda
package manager.

crema also depends on several Python packages:

- `numba <http://numba.pydata.org/>`_
- `numpy <https://numpy.org/>`_
- `pandas <https://pandas.pydata.org/>`_
- `pyteomics <https://pyteomics.readthedocs.io/en/latest/>`_




We recommend using `pip` to install crema. Missing dependencies will also
be installed automatically:

.. code-block:: bash

   $ pip3 install crema-ms

Basic Usage
-----------
Use **crema** from the Command Line
###################################

If your input files are in the mzTab or Crux tab-delimited formats, then simple crema analyses can be performed
straight from the command line!

Suppose your mzTab file is located at the directory "data/psms.mztab". Simply run the following command:

.. code-block:: bash

   $ crema data/psms.mztab

Suppose your Crux files are located at the directory "data/target_psms.txt" and "data/decoy_psms.txt".
Simply run the following command:

.. code-block:: bash

   $ crema data/target_psms.txt data/decoy_psms.txt

That's it. The software will then run the target-decoy competition FDR estimation method using information from these columns
to calculate confidence estimates for the given data.

Your results will be saved in your working directory as .txt files named "crema.psms.txt" and "crema.peptides.txt".
These files will contain an additional column ("crema q-value") that is appended to several columns
(specifically those that identify the psm, peptide sequence, and score) parsed from the input file.

For a full list of parameters, see the :doc:`Command Line Interface <cli>`.

Use **crema** as a Python package
###################################

Here is a simple demonstration of how to use crema as an API:

.. code-block:: Python

    >>> import crema
    >>> input_files = ["data/target_psms.txt", "decoy_psms/decoys.txt"]
    >>> psms = crema.read_crux(input_files)
    >>> results =  psms.assign_confidence(score_column="combined p-value", desc=True, eval_fdr=0.01, method="tdc")
    >>> results.to_txt(ouput_dir="example_output_dir", file_root="test", sep="\t", decoys=False)

Let's break this down and see what's really happening.

First, start up the Python interpreter:

.. code-block:: bash

   $ python3

Next, import crema as a package:

.. code-block:: Python

   >>> import crema

Call the :doc:`read_crux() <api/functions>` method and pass in the desired input files. Note that
the files "data/target_psms.txt" and "data/decoy_psms.txt" are already in the required Crux file
format. The :doc:`read_crux() <api/functions>` method will return a :doc:`dataset <api/dataset>` object
that we will save as "psms" in this example:

.. code-block:: Python

   >>> input_files = ["data/target_psms.txt", "decoy_psms/decoys.txt"]
   >>> psms = crema.read_crux(input_files)

Execute the desired FDR estimation method by calling the :doc:`assign_confidence <api/functions>` method on
the dataset object that we created above.
This operation will return a :doc:`confidence <api/confidence>` object that we will save as "results":

.. code-block:: Python

   >>> results =  psms.assign_confidence(score_column="combined p-value", desc=True, eval_fdr=0.01, method="tdc")

Note: The parameters passed here are optional and are only specified here for
demonstration. Further details can be found in the documentation for the :doc:`dataset <api/dataset>` class.

Confidence objects contain a :doc:`to_txt() <api/confidence>` method that allows you to write your results to a txt file.
Your results will be saved in your working directory (unless otherwise specified)
as .txt files named "crema.psms.txt" and "crema.peptides.txt".
These files will contain an additional column ("crema q-value") that is appended to several columns
(specifically those that identify the psm, peptide sequence, and score) parsed from the input file.

.. code-block:: Python

   >>> results.to_txt(ouput_dir="example_output_dir", file_root=None, sep="\t", decoys=False)

Note: The parameters passed here are optional and are only specified here for
demonstration. Further details can be found in the documentation for the :doc:`confidence <api/confidence>` class.

That's all there is to it! You have successfully used crema as an API to
calculate confidence estimates for your data.