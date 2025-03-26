import json
import os

import numpy as np
import pandas as pd
from Bio import SeqIO


class Config:
    def __init__(self, settings_filename: str):
        self.settings_filename = settings_filename
        self.settings = {}

    def _findfile(self, target) -> str | None:
        for dirpath, _, filename in os.walk(os.getcwd()):
            if target in filename:
                return os.path.join(dirpath, target)
        return None

    def _read_settings(self) -> str:
        settings_path = self._findfile(self.settings_filename)
        if not settings_path:
            return f"We couldn't find the file: {self.settings_filename}"

        with open(settings_path, "r") as s_f:
            data = json.load(s_f)
            if not data:
                return "The given file is empty"
            self.settings.update(**data)

        return ""

    def _read_score_matrix(self) -> str:
        score_matrix_filename = self.settings.get("score_matrix_file")
        if not score_matrix_filename:
            return f"'score_matrix' is missing from {self.settings_filename}"

        score_matrix_path = self._findfile(score_matrix_filename)
        if not score_matrix_path:
            return "We couldn't find the score matrix file: {score_matrix_filename}"

        df = pd.read_csv(score_matrix_path)
        # TODO: add check_score_matrix function
        if df.empty:
            return f"{score_matrix_filename} is empty"
        self.settings.update(
            {
                "score_matrix": np.array(df.values),
                "base_decode": {
                    char.strip().upper(): i for i, char in enumerate(df.columns)
                },
            }
        )
        return ""

    def _read_sequences(self, sequences_filename: str) -> str:
        if ".fasta" not in sequences_filename.lower():
            return "Sequences should be in FASTA format"

        sequences_path = self._findfile(sequences_filename)
        if not sequences_path:
            return f"We couldn't find the sequences file: {sequences_filename}"

        sequences = {
            record.id: str(record.seq)
            for record in SeqIO.parse(sequences_path, "fasta")
        }
        if len(sequences) < 2:
            return "There should be at least 2 sequences to align"

        self.settings.update({"sequences": sequences})
        return ""

    def parse_config(self, sequences_filename: str) -> str:
        error = self._read_settings()
        if error:
            return error

        error = self._read_score_matrix()
        if error:
            return error

        error = self._read_sequences(sequences_filename)
        return error
