"""
The :py:class:`PsmDataset` class is used to define a collection of peptide-spectrum matches
"""


class PsmDataset:
    """
    Store a collection of PSMs and their features.

    Parameters
        ----------
        data : pandas.DataFrame
            dataframe of PSMs with columns indicating spectrum, score, and target
        spectrum_col : str
            name of the column that identifies the psm
        score_col : str
            name of the column that defines the scores (p-values) of the psms
        target_col : str
            name of the column that indicates if a psm is a target/decoy


    Attributes
        ----------
        data : pandas.DataFrame
    """

    def __init__(self, data, spectrum_col, score_col, target_col):
        """
        Initialize a PsmDataset object.
        """

        self._data = data
        self.spectrum_col = spectrum_col
        self.score_col = score_col
        self.target_col = target_col

    @property
    def data(self):
        """The full data collection of PSMs as a :py:class:`pandas.DataFrame`."""
        return self._data

    def get_col(self, col_name):
        """The column specified by the col_name as a  :py:class:`pandas.DataFrame`."""
        return self._data.loc[:, col_name]
