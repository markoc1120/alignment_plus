import argparse
import sys

from config.config import Config
from scripts.msa_global import MSAAlignment
from scripts.pairwise_global import PairwiseAlignment
from scripts.utils import save_result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pairwise and MSA alignment with linear or affine gap costs.",
    )
    parser.add_argument(
        "-i", "--input", type=str, required=True, help="Your input fasta file"
    )
    parser.add_argument(
        "-at",
        "--alignment-type",
        choices=["pair", "multi"],
        required=True,
        help="Alignment type which can either be 'pair' for Pairwise or 'multi' for MSA",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        help="Output path to write the result",
    )
    parser.add_argument(
        "-gm",
        "--gap-model",
        choices=["linear", "affine"],
        default="linear",
        help="Gap penalty model to use, cost parameters can be set via settings JSON file (default: linear)",
    )
    parser.add_argument(
        "-s",
        "--settings",
        type=str,
        help="Name of the settings file (format should be JSON)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if args.settings:
        config = Config(args.settings)
        error = config.parse_config(args.input)
        if error:
            print(f"Error: {error}")

    # try:
    if args.alignment_type == "pair":
        model = PairwiseAlignment(config.settings)
    else:
        model = MSAAlignment(config.settings)
    error = model.check_initialization()
    score_matrix, alignments, err = model.align(args.gap_model.lower())
    if err:
        print(f"Error: {error}")
        sys.exit(1)

    if args.output:
        save_result(alignments, args.output)
        print(f"Aligned score: {score_matrix[-1, -1]}.")
    else:
        print(score_matrix)
        print(alignments)
    if error:
        print(f"Error: {error}")
    # except Exception as e:
    #     print(f"Error: {e}", file=sys.stderr)
    #     sys.exit(1)


if __name__ == "__main__":
    main()
