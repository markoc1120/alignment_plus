import os

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .msa_global import MSAAlignment
from .pairwise_global import PairwiseAlignment


def save_result(results: list[tuple[str, str]], output_path: str) -> None:
    if "/" in output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

    records = []
    for idx, (seq1, seq2) in enumerate(results, 1):
        record1 = SeqRecord(
            Seq(seq1), id=f"alignment_{idx}_seq1", description="First sequence"
        )
        record2 = SeqRecord(
            Seq(seq2), id=f"alignment_{idx}_seq2", description="Second sequence"
        )
        records.extend([record1, record2])
    SeqIO.write(records, output_path, "fasta")


def get_alignment_score(
    matrix: np.ndarray, alignment_model: PairwiseAlignment | MSAAlignment
):
    if isinstance(alignment_model, PairwiseAlignment):
        return matrix[-1, -1]
    elif isinstance(alignment_model, MSAAlignment):
        dims = matrix.shape
        score_idx = tuple(dim - 1 for dim in dims)
        return matrix[score_idx]
    else:
        raise ValueError("Unknown alignment model type")
