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
   vignettes/index.rst
   cli.rst
   api/index.rst
   notes.rst


Getting Started
---------------
**crema** produces confidence estimates for peptide detection in mass spectrometry proteomics experiments.
It takes files holding data regarding peptide-spectrum matches (PSMs) as input, executes the
desired estimation method, and produces confidence estimates of the PSMs as output.

Introduction
------------
One of the fundamental tasks in proteomics is detecting peptides from mass spectra that we identify through
Mass Spectrometry. We have tools that assign peptides to spectra, but unfortunately this matching is not 100% accurate -
meaning there is uncertainty about whether Peptide Spectrum Matches are real or false positives. We want to be able to
quantify this uncertainty so that we can be confident in our data and what it represents because this will further
ensure that expensive proteomics experiments use relevant and accurate data.

**crema** is a Python package that implements various methods to estimate false discovery rates (FDR) of peptide
detection in mass spectrometry proteomics experiments. Although there are many ways to estimate FDR, crema focuses on
methods that rely on the concept of target decoy competition. The sole purposes of crema is to do this, and to do this
well. As a result, we developed crema to be lightweight and flexible. It has very minimal dependencies and supports a
wide range of input and output formats. On top of that, it is extremely simple to use.

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

- `pandas <https://pandas.pydata.org/>`_


We recommend using `pip` to install crema. Missing dependencies will also
be installed automatically:

.. code-block:: bash

   $ pip3 install crema-ms

Basic Usage
-----------
Use **crema** from the Command Line
###################################

Simple crema analyses can be performed from the command line:

.. code-block:: bash

   $ crema data/single_basic.csv

That's it. Giving crema nothing but the input file will force it to search for
three specific column names: "combined p-value", "scan", and "target/decoy". It will then run the
Target-Decoy Competition FDR method using the information from these columns
to calculate confidence estimates for the given data.

Your results will be saved in your working directory as a
csv file named `crema.psm_results.txt`. This file will contain two additional columns
(False Discovery Rate and Q-Value) that are
appended to the initial few columns specified from the input file.

For a full list of parameters, see the :doc:`Command Line Interface <cli>`.

Use **crema** as a Python package
###################################

Here's a simple demonstration of how to use crema as an API:

.. code-block:: Python

   >>> import crema
   >>> psms = crema.read_file(["data/multi_target.csv", "data/multi_decoy.csv"])
   >>> results = crema.calculate_tdc(psms)
   >>> results.write_csv("save_to_here.txt")

Let's break this down and see what's really happening!


First, start up the Python interpreter:

.. code-block:: bash

   $ python3

Next, import crema as a package:

.. code-block:: Python

   >>> import crema

Call the read_file method and pass in the desired input files. In this example,
the files "data/multi_target.csv" and "data/multi_decoy.csv" are already in crux
format. Thus we do not need to specify non-default column names.
This will return a dataset object that we will save as "psms" in this example:

.. code-block:: Python

   >>> psms = crema.read_file(["data/multi_target.csv", "data/multi_decoy.csv"])

Execute the desired FDR estimation method by calling the "calculate_[algorithm]" method and
passing in the dataset object that we created above. This will return a result object that
we will save as "results" in this example:

.. code-block:: Python

   >>> results = crema.calculate_tdc(psms)

Result objects contain a "write_csv" method that allows you to write your result to a csv file.
Your results will be saved in your working directory (unless otherwise specified) as a
csv file named by the parameter you pass when calling the method.
This file will contain two additional columns
(False Discovery Rate and Q-Value) that are
appended to the initial few columns specified from the input file.

.. code-block:: Python

   >>> results.write_csv("save_to_here.txt")

That's all there is to it! You've successfully used crema as an API to
calculate confidence estimates for your data!