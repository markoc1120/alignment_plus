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
    for idx, sequences in enumerate(results, 1):
        for seq_num, sequence in enumerate(sequences, 1):
            record = SeqRecord(
                Seq(sequence),
                id=f"alignment_{idx}_seq{seq_num}",
                description=f"Sequence {seq_num}",
            )
            records.append(record)
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
