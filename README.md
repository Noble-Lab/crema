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
    >>> psms = crema.read_crux(input_files)
    >>> results =  psms.assign_confidence()
    >>> results.to_txt(ouput_dir="example_output_dir")
```

Check out our [documentation](hhttps://crema-ms.readthedocs.io) for more details
on how to make full use of crema.