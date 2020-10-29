.. crema documentation master file, created by
   sphinx-quickstart on Fri Sep 11 15:54:39 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to crema's documentation!
=================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


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

   $ pip3 install crema

Basic Usage
-----------
Use **crema** as a Python package
###################################

First, start up the Python interpreter:

.. code-block:: bash

   $ python3

Then calculate confidence estimate using crema:

.. code-block:: Python

   >>> import crema
   >>> psms = cream.read_file(["data/single.csv"], "scan", "p-value", "target")
   >>> result = crema.calculate_tdc(psms)
   >>> results.write_csv("save_to_here.csv")