from .alignment_base import AlignmentBase


class MSAAlignment(AlignmentBase):
    def __init__(self, settings):
        super().__init__(settings)

    def check_initialization(self) -> str:
        if len(self.sequences) > 2:
            return "Input FASTA file should have more than 2 sequences"
        return ""
