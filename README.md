<img src="https://raw.githubusercontent.com/Noble-Lab/crema/master/static/crema_logo.svg" width=300>
 
---

Confidence Estimation for Mass Spectrometry Proteomics

**crema** is a Python package that implements various methods to estimate false discovery rates (FDR)
in mass spectrometry proteomics experiments. crema focuses on
methods that rely on the concept of "target-decoy competition." The sole purpose of crema is to do decoy-based FDR
estimation, and to do it well. As a result, crema is lightweight and flexible. It has minimal dependencies and
supports a wide range of input and output formats. On top of that, it is extremely simple to use.

For more information, check out our
[documentation](https://crema-ms.readthedocs.io).  

## Installation  

crema requires Python 3.6+ and can be installed with pip:  

```
$ pip3 install crema-ms
```

## Basic Usage  

Before using crema, you need one or more files, each containing a collection of
peptide-spectrum matches (PSMs) in tab-delimited format. Note that crema defaults
to reading files via [crux](http://crux.ms/index.html) format, but can easily be
manipulated to accept files in formats that use differing column headers.

Simple crema calculations can be performed at the command line:

```Bash
$ crema data/tide-search.target.psms.txt data/tide-search.decoy.psms.txt
```

Alternatively, the Python API can be used to calculate confidence estimates in the Python
interpreter and affords greater flexibility:

```Python
    >>> import crema
    >>> input_files = ["data/tide-search.target.psms.txt", "data/tide-search.decoy.psms.txt"]
    >>> psms = crema.read_tide(input_files)
    >>> results =  psms.assign_confidence()
    >>> results.to_txt(ouput_dir="example_output_dir")
```

Check out our [documentation](hhttps://crema-ms.readthedocs.io) for more details
on how to make full use of crema.

## Disclaimer

This material was prepared as an account of work sponsored by an agency of the
United States Government.  Neither the United States Government nor the United
States Department of Energy, nor Battelle, nor any of their employees, nor any
jurisdiction or organization that has cooperated in the development of these
materials, makes any warranty, express or implied, or assumes any legal
liability or responsibility for the accuracy, completeness, or usefulness or any
information, apparatus, product, software, or process disclosed, or represents
that its use would not infringe privately owned rights.
Reference herein to any specific commercial product, process, or service by
trade name, trademark, manufacturer, or otherwise does not necessarily
constitute or imply its endorsement, recommendation, or favoring by the United
States Government or any agency thereof, or Battelle Memorial Institute. The
views and opinions of authors expressed herein do not necessarily state or
reflect those of the United States Government or any agency thereof.

<p align=center>PACIFIC NORTHWEST NATIONAL LABORATORY</p>
<p align=center><i>operated by</i></p>
<p align=center>BATTELLE</p>
<p align=center><i>for the</i></p>
<p align=center>UNITED STATES DEPARTMENT OF ENERGY</p>
<p align=center><i>under Contract DE-AC05-76RL01830</i></p>
