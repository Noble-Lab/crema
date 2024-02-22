"""The :py:class:`PsmDataset` class is used to define a collection of
peptide-spectrum matches.
"""

import logging

from .confidence import TdcConfidence
from .confidence import MixmaxConfidence
from .qvalues import tdc
from .utils import listify

LOGGER = logging.getLogger(__name__)


class PsmDataset:
    """Store a collection of peptide-spectrum matches (PSMs).

    Parameters
    ----------
    psms : pandas.DataFrame
        A :py:class:`pandas.DataFrame` of PSMs.
    target_column : str
        The column that indicates whether a PSM is a target or a decoy. This
        column should be boolean, where :code:`True` indicates a target and
        :code:`False` indicates a decoy.
    spectrum_columns : str or tuple of str
        One or more columns that together define a unique mass spectrum.
    score_columns : str or tuple of str, optional
        One or more columns that indicate scores by which crema can rank PSMs.
    peptide_column : str
        The column that defines a unique peptide. Modifications should be
        indicated either in square brackets :code:`[]` or parentheses
        :code:`()`. The exact modification format within these entities does
        not matter, so long as it is consistent.
    protein_columns : str
        The column that defines a unique protein.
    protein_delim : str
        The string delimiter that is needed to separate multiple proteins found
        in the protein column.
    peptide_pairing: dict[str, str]
        A map of target and decoy peptide sequence pairings to be used for TDC.
        This should be in the form {key=target_sequence:value=decoy_sequence}
        where decoy sequences are shuffled versions of target sequences.
    copy_data : bool, optional
        If true, a deep copy of the data is created. This uses more memory, but
        is safer because it prevents accidental modification of the underlying
        data. This argument only has an effect when `pin_files` is a
        :py:class:`pandas.DataFrame`

    Attributes
    ----------
    data : pandas.DataFrame
    spectrum_columns : list of str
    score_columns : list of str
    target_column : str
    peptide_column : str
    protein_column : str
    protein_delim : str
    methods : dict
    peptide_pairing : dict
    """

    methods = {"tdc": TdcConfidence, "mixmax": MixmaxConfidence}

    def __init__(
        self,
        psms,
        target_column,
        spectrum_columns,
        score_columns,
        peptide_column,
        protein_column,
        protein_delim,
        peptide_pairing=None,
        copy_data=True,
    ):
        """Initialize a PsmDataset object."""
        self.score_columns = listify(score_columns)
        self._spectrum_columns = listify(spectrum_columns)
        self._target_column = target_column
        self._peptide_column = peptide_column
        self._protein_column = protein_column
        self._protein_delim = protein_delim
        self._peptide_pairing = peptide_pairing

        fields = sum(
            [
                self._spectrum_columns,
                self.score_columns,
                [self._target_column],
                [self._peptide_column],
                [self._protein_column],
            ],
            [],
        )
        self._data = psms.copy(deep=copy_data).loc[:, fields]
        self._data[target_column] = self._data[target_column].astype(bool)
        self._num_targets = self.targets.sum()
        self._num_decoys = (~self.targets).sum()

        if self.data.empty:
            raise ValueError("No PSMs were detected.")

        if not self._num_decoys:
            raise ValueError("No decoy PSMs were detected.")

        if not self._num_targets:
            raise ValueError("No target PSMs were detected.")

    @property
    def columns(self):
        """The columns of the PSM :py:class:`pandas.DataFrame`"""
        return self._data.columns

    @property
    def data(self):
        """The collection of PSMs as a :py:class:`pandas.DataFrame`."""
        return self._data.copy()

    @property
    def spectra(self):
        """The mass spectrum identifiers as a :py:class:`pandas.DataFrame`."""
        return self[self._spectrum_columns]

    @property
    def peptides(self):
        """The peptides as a :py:class:`pandas.Series`."""
        return self[self._peptide_column]

    @property
    def proteins(self):
        """The proteins as a :py:class:`pandas.Series`."""
        return self[self._protein_column]

    @property
    def protein_delim(self):
        """The delimiter to split protein IDs as a string."""
        return self[self._protein_delim]

    @property
    def scores(self):
        """The scores for each PSM as a :py:class:`pandas.DataFrame`."""
        return self[self.score_columns]

    @property
    def targets(self):
        """An array indicating whether each PSM is a target"""
        return self[self._target_column].values

    @property
    def peptide_pairing(self):
        """A dictionary containing target/decoy peptide pairs"""
        return self._peptide_pairing

    def __getitem__(self, column):
        """Return the specified column"""
        return self._data.loc[:, column]

    def assign_confidence(
        self,
        score_column=None,
        threshold=0.01,
        pep_fdr_type="psm-peptide",
        prot_fdr_type="best",
        desc=None,
        eval_fdr=0.01,
        method="tdc",
    ):
        """Assign confidence estimates to this collection of peptide-spectrum matches.

        Parameters
        ----------
        score_column : str, optional
            The score by which to rank the PSMs for confidence estimation. If
            :code:`None`, the score that yields the most PSMs at the specified
            false discovery rate threshold (`eval_fdr`), will be used.
        threshold : float or "q-value", optional
            The FDR threshold for accepting discoveries. Default is 0.01. If
            "q-value" is chosen, then "accept" column is replaced with
            "crema q-value".
        pep_fdr_type : {"psm-only","peptide-only",psm-peptide"}, optional
            The method for crema to use when calculating peptide level confidence
            estimates. Default is "psm-peptide".
        prot_fdr_type : {"best", "combine"}, optional
            The method for crema to use when calculating protein level confidence
            estimates. Default is "best".
        desc : bool, optional
            True if higher scores better, False if lower scores are better.
            If None, crema will try both and use the
            choice that yields the most PSMs at the specified false discovery
            rate threshold (`eval_fdr`). If `score_column` is :code:`None`,
            this parameter is ignored.
        eval_fdr : float, optional
            The false discovery rate threshold used to evaluate the best
            `score_column` and `desc` to choose. This should range from 0 to 1.
        method : {"tdc"}, optional
            The method for crema to use when calculating the confidence estimates.

        Returns
        -------
        Confidence object
            The confidence estimates for this PsmDataset.
        """
        if score_column is None:
            score_column, _, desc = self.find_best_score(eval_fdr)

        conf = self.methods[method](
            psms=self,
            score_column=score_column,
            desc=desc,
            eval_fdr=eval_fdr,
            pep_fdr_type=pep_fdr_type,
            prot_fdr_type=prot_fdr_type,
            threshold=threshold,
        )

        return conf

    def find_best_score(self, eval_fdr=0.01):
        """Find the best score for this collection of PSMs

        Try each of the available score columns, determining how many PSMs
        are detected below the provided false discovery rate threshold. The
        best score is the one that returns the most.

        Parameters
        ----------
        eval fdr : float
            The false discovery rate threshold used to find the best score.

        Returns
        -------
        best_score : str
            The best score.
        n_passing : int
            The number of PSMs that meet the specified FDR threshold.
        desc : bool
            True if higher scores better, False if lower scores are better.
        """
        best_score = None
        best_passing = 0
        for desc in (True, False):
            qvals = self.scores.apply(tdc, target=self.targets, desc=desc)
            num_passing = (qvals <= eval_fdr).sum()
            feat_idx = num_passing.idxmax()
            num_passing = num_passing[feat_idx]
            if num_passing > best_passing:
                best_passing = num_passing
                best_score = feat_idx
                best_desc = desc

        if best_score is None:
            raise RuntimeError("No PSMs were found below the 'eval_fdr'.")

        return best_score, best_passing, best_desc

    def set_protein_column(self, new_protein_column):
        """Replaces current protein column with input protein column

        Parameters
        ----------
        new_protein_column : pandas.Series

        Returns
        -------

        """
        self._data[self._protein_column] = new_protein_column
        return

    def set_peptide_column(self, new_peptide_column):
        """Replaces current peptide column with input peptide column

        Parameters
        ----------
        new_peptide_column : pandas.Series

        Returns
        -------

        """
        self._data[self._peptide_column] = new_peptide_column
        return
