Python API
==========
The Python API enables maximum flexibility using crema.

Read PSMs easily using the :py:func:`~mokapot.read_crux()` or
:py:func:`~mokapot.read_mztab()` functions for files in the Crux
tab-delimited format or MzTab format, respectively. Alternatively,
more generic tab delimited files can be read using the :py:func:`~mokapot.read_txt()`.

Once read, a collection of PSMs is stored in a :py:class:`~crema.dataset.PsmDataset` object.
Calling the :py:func:`~crema.assign_confidence()` function will calculate confidence estimates
for the respective psm data using the specified confidence estimate method.

This will produce a :py:class:`~crema.Confidence` object that stores the calculated
confidence estimates via False Discovery Rates and Q-Values. Results stored in a
confidence object can be exported by calling the :py:func:`~crema.to_txt()` function.


.. toctree::
    :maxdepth: 1
    :hidden:
    :titlesonly:

    Overview <self>
    functions.rst
    dataset.rst
    confidence.rst


Functions
---------
.. currentmodule:: crema

Primary Functions
*****************
.. autosummary::
    :nosignatures:

    assign_confidence

Parsers
*****************
.. autosummary::
    :nosignatures:

    read_crux
    read_mztab
    read_txt

Writers
*****************
.. autosummary::
    :nosignatures:

    to_txt

Dataset
---------
.. currentmodule:: crema.dataset
.. autosummary::
    :nosignatures:

    PsmDataset

Confidence
---------
.. currentmodule:: crema.confidence
.. autosummary::
    :nosignatures:

    Confidence


