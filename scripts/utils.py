import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
