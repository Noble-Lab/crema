<img src="https://github.com/Noble-Lab/crema/blob/master/docs/_static/crema_logo_caramel_light.png" width=300>  

---

Confidence Estimation for Mass Spectrometry Proteomics

**crema** is a Python package that implements various methods to estimate false discovery rates (FDR) of peptide
detection in mass spectrometry proteomics experiments. Although there are many ways to estimate FDR, crema focuses on
methods that rely on the concept of target decoy competition. The sole purposes of crema is to do this, and to do this
well. As a result, we developed crema to be lightweight and flexible. It has very minimal dependencies and supports a
wide range of input and output formats. On top of that, it is extremely simple to use.

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
$ crema data/single_basic.csv
```

Alternatively, the Python API can be used to calculate confidence estimates in the Python
interpreter and affords greater flexibility:

```Python
    >>> import crema-ms
    >>> psms = crema.read_file(["data/multi_target.csv", "data/multi_decoy.csv"])
    >>> results = crema.calculate_tdc(psms)
    >>> results.write_csv("save_to_here.txt")
```

Check out our [documentation](hhttps://crema-ms.readthedocs.io) for more details
on how to make full use of crema.