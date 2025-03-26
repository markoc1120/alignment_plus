import numpy as np

from config.constants import AFFINE_GAP_MODEL, LINEAR_GAP_MODEL

from .alignment_base import AlignmentBase


class PairwiseAlignment(AlignmentBase):
    def __init__(self, settings):
        super().__init__(settings)

    def _init_C(self, seq1_length: int, seq2_length: int) -> np.ndarray:
        C = np.zeros((seq1_length, seq2_length))
        C[0, :] = range(seq2_length)
        C[:, 0] = range(seq1_length)
        C[0, :] = C[0, :] * self.gapopen
        C[:, 0] = C[:, 0] * self.gapopen
        return C

    def _calc_cost_linear(self, seq1: str, seq2: str) -> np.ndarray:
        # init C table based on sequences lengths
        seq1_length, seq2_length = len(seq1) + 1, len(seq2) + 1
        C = self._init_C(seq1_length, seq2_length)

        # fill C row by row
        for s1_idx in range(1, seq1_length):
            for s2_idx in range(1, seq2_length):
                match_cost = self.score_matrix[
                    self.base_decode[seq1[s1_idx - 1].upper()],
                    self.base_decode[seq2[s2_idx - 1].upper()],
                ]
                # chosse minimum cost for the three cases
                C[s1_idx, s2_idx] = min(
                    C[s1_idx - 1, s2_idx] + self.gapopen,  # deletion
                    C[s1_idx, s2_idx - 1] + self.gapopen,  # insertion
                    C[s1_idx - 1, s2_idx - 1] + match_cost,  # match/mismatch
                )
        return C

    def _backtrack_linear(
        self, C: np.ndarray, seq1: str, seq2: str
    ) -> list[tuple[str, str]]:
        s1_idx, s2_idx = len(seq1), len(seq2)
        align1 = align2 = ""
        results = []

        def trace_alignments(s1_idx, s2_idx, align1, align2):
            # check if we have depleted our indexes, save the two alignments
            if s1_idx == s2_idx == 0:
                results.append((align1, align2))
                return

            # check if we found a match/mismatch -> add recursive call to stack
            if (
                s1_idx > 0
                and s2_idx > 0
                and C[s1_idx, s2_idx]
                == C[s1_idx - 1, s2_idx - 1]
                + self.score_matrix[
                    self.base_decode[seq1[s1_idx - 1].upper()],
                    self.base_decode[seq2[s2_idx - 1].upper()],
                ]
            ):
                trace_alignments(
                    s1_idx - 1,
                    s2_idx - 1,
                    seq1[s1_idx - 1] + align1,
                    seq2[s2_idx - 1] + align2,
                )
            # check if we found a deletion -> add recursive call to stack
            if s1_idx > 0 and C[s1_idx - 1, s2_idx] + self.gapopen == C[s1_idx, s2_idx]:
                trace_alignments(
                    s1_idx - 1, s2_idx, seq1[s1_idx - 1] + align1, "-" + align2
                )
            # check if we found an insertion -> add recursive call to stack
            if s2_idx > 0 and C[s1_idx, s2_idx - 1] + self.gapopen == C[s1_idx, s2_idx]:
                trace_alignments(
                    s1_idx, s2_idx - 1, "-" + align1, seq2[s2_idx - 1] + align2
                )

        trace_alignments(s1_idx, s2_idx, align1, align2)
        return results

    def _calc_cost_affine(self, seq1: str, seq2: str) -> np.ndarray:
        seq1_length, seq2_length = len(seq1) + 1, len(seq2) + 1

        S = np.full((seq1_length, seq2_length), float("inf"))
        D = np.full((seq1_length, seq2_length), float("inf"))
        I = np.full((seq1_length, seq2_length), float("inf"))
        S[0, 0] = 0

        # fill S, D, I row by row
        for i in range(seq1_length):
            for j in range(seq2_length):
                if i == j == 0:
                    continue

                if i > 0 and j > 0:
                    match_score = self.score_matrix[
                        self.base_decode[seq1[i - 1].upper()],
                        self.base_decode[seq2[j - 1].upper()],
                    ]
                    s_match = S[i - 1, j - 1] + match_score

                # calc D(i,j)
                if i > 0:
                    d_scores = []
                    if j >= 0:
                        d_scores.extend(
                            [
                                S[i - 1, j] + self.gapopen + self.gapextend,
                                D[i - 1, j] + self.gapextend,
                            ]
                        )
                    if d_scores:
                        D[i, j] = min(d_scores)

                # calc I(i,j)
                if j > 0:
                    i_scores = []
                    if i >= 0:
                        i_scores.extend(
                            [
                                S[i, j - 1] + self.gapopen + self.gapextend,
                                I[i, j - 1] + self.gapextend,
                            ]
                        )
                    if i_scores:
                        I[i, j] = min(i_scores)

                # calc S(i,j)
                s_scores = []
                if i > 0 and j > 0:
                    s_scores.append(s_match)
                if i > 0 and j >= 0:
                    s_scores.append(D[i, j])
                if i >= 0 and j > 0:
                    s_scores.append(I[i, j])
                if s_scores:
                    S[i, j] = min(s_scores)
        return S

    def _backtrack_affine(
        self, S: np.ndarray, seq1: str, seq2: str
    ) -> list[tuple[str, str]]:
        s1_idx, s2_idx = len(seq1), len(seq2)
        align1 = align2 = ""
        results = []

        def trace_alignments(s1_idx, s2_idx, align1, align2):
            # check if we have depleted our indexes, save the two alignments
            if s1_idx == s2_idx == 0:
                results.append((align1, align2))
                return

            curr_score = S[s1_idx, s2_idx]
            # check if we found a match/mismatch -> add recursive call to stack
            if (
                s1_idx > 0
                and s2_idx > 0
                and curr_score
                == S[s1_idx - 1, s2_idx - 1]
                + self.score_matrix[
                    self.base_decode[seq1[s1_idx - 1].upper()],
                    self.base_decode[seq2[s2_idx - 1].upper()],
                ]
            ):
                trace_alignments(
                    s1_idx - 1,
                    s2_idx - 1,
                    seq1[s1_idx - 1] + align1,
                    seq2[s2_idx - 1] + align2,
                )

            # try deletion(s) and insertion(s)
            k = 1
            while s1_idx >= k or s2_idx >= k:
                # check if we found deletion(s) -> add recursive call to stack
                if s1_idx >= k and curr_score == S[
                    s1_idx - k, s2_idx
                ] + self.gapopen + (k * self.gapextend):
                    trace_alignments(
                        s1_idx - k,
                        s2_idx,
                        seq1[s1_idx - k : s1_idx] + align1,
                        "-" * k + align2,
                    )

                # check if we found insertion(s) -> add recursive call to stack
                if s2_idx >= k and curr_score == S[
                    s1_idx, s2_idx - k
                ] + self.gapopen + (k * self.gapextend):
                    trace_alignments(
                        s1_idx,
                        s2_idx - k,
                        "-" * k + align1,
                        seq2[s2_idx - k : s2_idx] + align2,
                    )

                k += 1

        trace_alignments(s1_idx, s2_idx, align1, align2)
        return results

    def check_initialization(self) -> str:
        if len(self.sequences) != 2:
            return "Input FASTA file should exactly have 2 sequences"
        return ""

    def pairwise_alignment_linear(
        self, seq1: str, seq2: str
    ) -> tuple[np.ndarray, list[tuple[str, str]], str]:
        C = self._calc_cost_linear(seq1, seq2)
        alignments = self._backtrack_linear(C, seq1, seq2)
        return C, alignments, ""

    def pairwise_alignment_affine(
        self, seq1: str, seq2: str
    ) -> tuple[np.ndarray, list[tuple[str, str]], str]:
        S = self._calc_cost_affine(seq1, seq2)
        alignments = self._backtrack_affine(S, seq1, seq2)
        return S, alignments, ""

    def align(
        self, gap_model: str
    ) -> tuple[np.ndarray | None, list[tuple[str, str]] | None, str]:
        seq1, seq2 = self.sequences.values()
        if gap_model == LINEAR_GAP_MODEL:
            return self.pairwise_alignment_linear(seq1, seq2)
        if gap_model == AFFINE_GAP_MODEL:
            return self.pairwise_alignment_affine(seq1, seq2)
        return None, None, ""
