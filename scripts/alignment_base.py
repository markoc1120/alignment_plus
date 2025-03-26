import numpy as np


class AlignmentBase:
    def __init__(self, settings):
        self.settings = settings

    @property
    def gapopen(self) -> int | None:
        return self.settings.get("gapopen")

    @property
    def gapextend(self) -> int | None:
        return self.settings.get("gapextend")

    @property
    def score_matrix(self) -> np.ndarray:
        return self.settings.get("score_matrix", {})

    @property
    def base_decode(self) -> dict[str, int]:
        return self.settings.get("base_decode", {})

    @property
    def sequences(self) -> dict[str, str]:
        return self.settings.get("sequences", {})
