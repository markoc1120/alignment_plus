import numpy as np

from config.constants import LINEAR_GAP_MODEL

from .alignment_base import AlignmentBase


class MSAAlignment(AlignmentBase):

    def __init__(self, settings):
        super().__init__(settings)

    def _init_score_matrixes(
        self, dimensions: tuple[int, ...], sequences: list[str]
    ) -> np.ndarray:
        score_matrixes = np.zeros(dimensions, dtype=int)
        for idx in np.ndindex(dimensions):
            _sum = sum(idx)
            if _sum == 0:
                continue
            score_matrixes[idx] = _sum * self.gapopen
        return score_matrixes

    def _calculate_score_for_column(
        self, prev_idx: tuple[int, ...], delta: tuple[int, ...], sequences: list[str]
    ) -> int:
        column_score = 0
        num_sequences = len(sequences)

        for i in range(num_sequences):
            for j in range(i + 1, num_sequences):
                if delta[i] == delta[j] == 1:
                    column_score += self.score_matrix[
                        self.base_decode[sequences[i][prev_idx[i]].upper()],
                        self.base_decode[sequences[j][prev_idx[j]].upper()],
                    ]
                else:
                    column_score += self.gapopen
        return column_score

    def _get_possible_scores_and_prev_indexes(
        self, idx: tuple[int, ...], score_matrixes: np.ndarray, sequences: list[str]
    ) -> list[tuple[tuple[int, ...], int]]:
        num_sequences = len(sequences)
        scores = []

        for delta in np.ndindex((2,) * num_sequences):
            if sum(delta) == 0:
                continue
            prev_idx = tuple(np.subtract(idx, delta))
            if any(idx < 0 for idx in prev_idx):
                continue

            column_score = self._calculate_score_for_column(prev_idx, delta, sequences)
            scores.append((prev_idx, score_matrixes[prev_idx] + column_score))
        return scores

    def _calculate_scores(self, score_matrixes: np.ndarray, sequences: list[str]):
        for idx in np.ndindex(score_matrixes.shape):
            if sum(idx) == 0:
                continue
            possible_scores = [
                score
                for _, score in self._get_possible_scores_and_prev_indexes(
                    idx, score_matrixes, sequences
                )
            ]
            score_matrixes[idx] = max(possible_scores)

    def _backtrack_linear(
        self, score_matrixes: np.ndarray, sequences: list[str]
    ) -> list[str]:
        aligned_sequences = ["" for _ in sequences]

        idx = tuple(len(seq) for seq in sequences)
        while any(i > 0 for i in idx):
            predecessors = self._get_possible_scores_and_prev_indexes(
                idx, score_matrixes, sequences
            )
            best_prev_idx, _ = min(predecessors, key=lambda x: x[1])
            delta = tuple(np.subtract(idx, best_prev_idx))
            for i, d in enumerate(delta):
                if d == 1:
                    aligned_sequences[i] = (
                        sequences[i][best_prev_idx[i]] + aligned_sequences[i]
                    )
                else:
                    aligned_sequences[i] = "-" + aligned_sequences[i]
            idx = best_prev_idx
        return aligned_sequences

    def check_initialization(self) -> str:
        if len(self.sequences) <= 2:
            return "Input FASTA file should have more than 2 sequences"
        return ""

    def msa_alignment_linear(
        self, sequences: list[str]
    ) -> tuple[np.ndarray, list[tuple[str, str]], str]:
        dimensions = tuple(len(seq) + 1 for seq in sequences)
        score_matrixes = self._init_score_matrixes(dimensions, sequences)
        self._calculate_scores(score_matrixes, sequences)
        alignments = self._backtrack_linear(score_matrixes, sequences)
        return score_matrixes, alignments, ""

    def align(self, gap_model: str) -> tuple[np.ndarray | None, list[str] | None, str]:
        sequences = list(self.sequences.values())
        if gap_model == LINEAR_GAP_MODEL:
            return self.msa_alignment_linear(sequences)
        return None, None, ""
