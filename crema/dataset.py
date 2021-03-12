"""
The :py:class:`PsmDataset` class is used to define a collection of peptide-spectrum matches
"""


class PsmDataset:
    """
    Store a collection of PSMs and their features.

    Parameters
        ----------
        data : pandas.DataFrame
            Dataframe of PSMs with columns indicating spectrum, score, and target
        sequence_col : str or tuple of str, optional
            One or more column names that identify the peptide sequence.
        spectrum_col : str or tuple of str, optional
            One or more column names that identify the psm.
        score_col : str or tuple of str, optional
            One or more column names that define the scores (p-values) of the psms.
        target_col : str
            Name of the column that indicates if a psm is a target/decoy.


    Attributes
        ----------
        data : pandas.DataFrame
    """

    def __init__(
        self,
        data,
        sequence_col,
        spectrum_col,
        score_col,
        target_col,
        protein_map=None,
    ):
        """
        Initialize a PsmDataset object.
        """

        self._data = data
        self.sequence_col = sequence_col
        self.spectrum_col = spectrum_col
        self.score_col = score_col
        self.target_col = target_col
        self.protein_map = protein_map
        self.results = {
            "psm": [],
            "peptide": [],
            "protein": [],
        }

    @property
    def data(self):
        """The full data collection of PSMs as a :py:class:`pandas.DataFrame`."""
        return self._data

    def get_col(self, col_name):
        """The column specified by the col_name as a  :py:class:`pandas.DataFrame`."""
        return self._data.loc[:, col_name]

    def write_results(self, output_dir=None, file_root=None):
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
        # for i in range(len(self.results)):
        #     self.results[i].write_file(output_dir, file_root=str(i))
        for key in self.results:
            for i in range(len(self.results.get(key))):
                self.results.get(key)[i].write_file(
                    output_dir=output_dir,
                    file_root=file_root,
                    indicator=str(i),
                )
