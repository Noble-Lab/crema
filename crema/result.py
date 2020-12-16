"""
The :py:class:`Result` class is used to define a collection of peptide-spectrum matches with calculated
False Discovery Rates and Q-Values.
"""

import os


class Result:
    """
    Store a collection of PSMs, their features, and their FDR/Q-Values

    Parameters
        ----------
        data : pandas.DataFrame
            dataframe of PSMs with columns indicating spectrum, score, and target
        spectrum_col : str or tuple of str
            one or more column names that identifies the psm
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

        self._data = data.reset_index(drop=True)
        self.spectrum_col = spectrum_col
        self.score_col = score_col
        self.target_col = target_col
        self.fdr_col = "FDR"
        self.q_val_col = "Q_Value"

    @property
    def data(self):
        """The full data collection of PSMs as a :py:class:`pandas.DataFrame`."""
        return self._data

    def get_col(self, col_name):
        """The column specified by the col_name as a  :py:class:`pandas.DataFrame`."""
        return self._data.loc[:, col_name]

    def write_file(self, output_dir=None, file_root=None):
        """
        Exports the data as a .txt file with the suffix "crema.psm_results.txt".

        Parameters
        ----------
        output_dir : str, optional
            The directory in which to save the files. Defaults to the current working directory if not specified.
        file_root : str, optional
            A prefix concatenated to the output result file. Defaults to none.

        Returns
        -------
        str
            The file path to the exported results file
        """
        out_file = "crema.psm_results.txt"
        if output_dir is None:
            output_dir = os.getcwd()
        if file_root is not None:
            out_file = file_root + out_file
        file_path = os.path.join(output_dir, out_file)
        self.data.to_csv(file_path)
        return file_path
