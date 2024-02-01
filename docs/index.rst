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
   faq.rst
   Changelog <CHANGELOG.md>


Getting Started
---------------
**Crema** produces confidence estimates for peptide detection in mass spectrometry proteomics experiments.
It takes as input one or more files containing peptide-spectrum matches (PSMs), executes the
desired estimation method, and produces as output confidence estimates at the PSM, peptide and protein levels.

Introduction
------------
One of the fundamental tasks in mass spectrometry proteomics is detecting peptides on the basis of the observed mass
spectra. Many tools exist to assign peptides to spectra, but unfortunately this matching is never 100% accurate,
meaning that there is always uncertainty about whether a given peptide detection is correct or a false positive. We want to be able to
quantify this uncertainty so that we can be confident in our conclusions and ensure that expensive
downstream validation experiments use relevant and accurate data.

Crema is a Python package that implements various methods to estimate false discovery rates (FDRs)
in mass spectrometry proteomics experiments. Crema focuses on
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

Crema also depends on several Python packages:

- `numba <http://numba.pydata.org/>`_
- `numpy <https://numpy.org/>`_
- `pandas <https://pandas.pydata.org/>`_
- `lxml <https://lxml.de/>`_
- `pyteomics <https://pyteomics.readthedocs.io/en/latest/>`_


We recommend using `pip` to install crema. Missing dependencies will also
be installed automatically:

.. code-block:: bash

   $ pip3 install crema-ms

Basic Usage
-----------
Use **crema** from the Command Line
###################################

If your input files are in one of crema's supported file formats, such mzTab or
Tide tab-delimited, then simple crema analyses can be performed
straight from the command line.

Suppose your mzTab file is located at the directory "data/psms.mztab". Simply run the following command:

.. code-block:: bash

   $ crema data/psms.mztab

Alternatively, if your Tide files are located in "data/target_psms.txt" and "data/decoy_psms.txt",
then you would run the following command:

.. code-block:: bash

   $ crema data/target_psms.txt data/decoy_psms.txt

That's it. The software will run the target-decoy competition FDR estimation method using information from your files
to calculate confidence estimates for the given data.

Your results will be saved in your working directory as .txt files named
"crema.psms.txt", "crema.peptides.txt", "crema.proteins.txt", and "crema.protein_groups.txt".
These files will contain an additional column ("crema q-value") that is appended to several columns
parsed from the input file. A log file will also be saved "crema.log.txt".

For a full list of parameters, see the :doc:`Command Line Interface <cli>`.

Use **crema** as a Python Package
###################################

Here is a simple demonstration of how to use crema as an API:

.. code-block:: Python

    >>> import crema
    >>> input_files = ["data/example_psms_target.txt",
    >>>                "decoy_psms/example_psms_decoy.txt"]
    >>> pairing_file = "pairing_file.txt"
    >>> psms = crema.read_tide(input_files, pairing_file_name=pairing_file)
    >>> results = psms.assign_confidence(score_column="combined p-value",
    >>>           pep_fdr_type="psm-peptide", threshold=0.01)
    >>> results.to_txt(output_dir="example_output_dir", file_root="test", sep="\t", decoys=False)

Let's break this down and see what's really happening.

First, start up the Python interpreter:

.. code-block:: bash

   $ python3

Next, import crema as a package:

.. code-block:: Python

   >>> import crema

Call the :doc:`read_tide() <api/functions>` method and pass in the desired input
files. The files "data/example_psms_target.txt" and
"data/example_psms_decoys.txt" contains PSMs
and are in the required Tide file format. In addition, the pairing_file is an
optional argument that explicitly pairs target and decoy peptides.
The :doc:`read_tide() <api/functions>` method will return a :doc:`dataset <api/dataset>` object
that we will save as "psms" in this example:

.. code-block:: Python

   >>> input_files = ["data/example_psms_target.txt",
   >>>                "decoy_psms/example_psms_target.txt"]
   >>> pairing_file = "pairing_file.txt"
   >>> psms = crema.read_tide(input_files. pairing_file_name=pairing_file)

Note that you can replace :doc:`read_tide() <api/functions>` with other methods
such as :doc:`read_txt() <api/functions>` and :doc:`read_msgf()
<api/functions>`. Also note that while, in this example, the target and decoy PSMs are
separate files, they can combined together and passed as a single file.

Execute the desired FDR estimation method by calling the :doc:`assign_confidence <api/functions>` method on
the dataset object that we created above.
This operation will return a :doc:`confidence <api/confidence>` object that we will save as "results":

.. code-block:: Python

   >>> results =  psms.assign_confidence(score_column="combined p-value",
   >>>            pep_fdr_type="psm-peptide", threshold=0.01)

Note that the parameters passed here are optional and are only specified here for
demonstration. In this command, the score_column argument specifies which column
contains the PSM scores and the threshold argument specifies the FDR threshold
for which to accept a PSM. A full list of the parameters and further details
can be found in the documentation for the :doc:`dataset <api/dataset>` class.

Also note that the pep_fdr_type argument denotes the method used to estimate peptide-level
FDR. This argument supports three options: psm-only, peptide-only, and psm-peptide. A pairing
file is required to run the peptide-only or psm-peptide options. Also, note that
peptide-only requires a separate target and decoy database search. If
peptide-only is used in conjunction with a concatenated target-decoy search, then
it becomes equivalent to psm-peptide.

Confidence objects contain a :doc:`to_txt() <api/confidence>` method that allows
you to write your results to a text file. Your results will be saved in your
working directory (unless otherwise specified) as text files named
"crema.psms.txt", "crema.peptides.txt", "crema.proteins.txt",
and "crema.protein_groups.txt". These files will contain an additional column
("crema q-value") that is appended to several columns parsed from the input file.

.. code-block:: Python

   >>> results.to_txt(output_dir="example_output_dir", file_root=None, sep="\t", decoys=False)

Note that the parameters passed here are optional and are only specified here for
demonstration. Further details can be found in the documentation for the :doc:`confidence <api/confidence>` class.

That's all there is to it. You have successfully used crema as an API to
calculate confidence estimates for your data.

Plotting performance curves
###################################
Instead of outputting a list of discoveries, Crema can also be used for
benchmarking. One common way to perform this task is to plot the number of
detections as a function of FDR threshold. We have opted not to create a plot
function within Crema as it would add a dependency and there are existing user
friendly libraries for plotting. In lieu of a plot function, we provide example
code to show users how to create such a plot.

Below is a simple demonstration of how to run crema and plot the number of
detections as a function of FDR threshold. 

However, we emphasize that this type of PSM-level plot should only be used
for benchmarking studies and should not be used for experiments that
result in some biological claim.

.. code-block:: Python

    >>> import crema
    >>> from matplotlib import pyplot as plt
    >>> import seaborn as sns
    >>>
    >>> input_files = ["data/example_psms_target.txt",
    >>>                "decoy_psms/example_psms_decoy.txt"]
    >>> pairing_file = "pairing_file.txt"
    >>> psms = crema.read_tide(input_files, pairing_file_name=pairing_file)
    >>> results = psms.assign_confidence(score_column="combined p-value",
    >>>           pep_fdr_type="psm-peptide", threshold="q-value")
    >>>
    >>> # extract peptide-level q-values
    >>> detections = results.confidence_estimates['peptides']['crema q-value']
    >>> # remove detections with a q-value > 0.11
    >>> detections = detections[detections <= 0.11]
    >>> 
    >>> # create plot
    >>> fig,ax = plt.subplots()
    >>> sns.ecdfplot(detections,ax=ax,stat='count') 
    >>> plt.xlim(0,.10)
    >>> plt.ylabel("# of confident IDs")
    >>> plt.xlabel("crema q-value")
    >>> plt.tight_layout()
    >>> plt.savefig("plot.pdf")
    >>> plt.close() 

The difference between this code snippet and the original code snippet above 
is that we have also imported the plotting packages matplotlib
and seaborn. In addition, we have given the threshold argument in the 
:doc:`assign_confidence <api/functions>` call a value of
"q-value" instead of a float. This is required for crema to output q-values.
Finally, we have added eight lines of code to create and save the performance
curve. 

Supported Database Search Engines
###################################

Crema currently supports output generated from Tide, MSGF+, MSAmanda, Comet, MSFragger.

In addition, crema supports input files from any search engine that are in the
following file formats: mzTab, pepXML, and generic tab-delimited text files.
