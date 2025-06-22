import argparse
import sys

from config.config import Config
from scripts.msa_global import MSAAlignment
from scripts.pairwise_global import PairwiseAlignment
from scripts.utils import get_alignment_score, save_result


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

    err = model.check_initialization()
    if err:
        print(f"Error: {err}")
        sys.exit(1)

    score_matrix, alignments, err = model.align(args.gap_model.lower())
    if err:
        print(f"Error: {error}")
        sys.exit(1)

    try:
        alignment_score = get_alignment_score(score_matrix, model)
    except Exception as e:
        print(f"Error calculating score: {e}")
        sys.exit(1)

    print(score_matrix)
    if args.output:
        save_result(alignments, args.output)
    else:
        print(alignments)
    print(f"Aligned score: {alignment_score}.")
    # except Exception as e:
    #     print(f"Error: {e}", file=sys.stderr)
    #     sys.exit(1)


if __name__ == "__main__":
    main()
