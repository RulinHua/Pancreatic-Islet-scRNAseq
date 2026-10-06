#!/usr/bin/env python3
"""Export pySCENIC regulon targets to a tab-separated table."""

import argparse
import csv
import os
from pyscenic.cli.utils import load_signatures


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--scenic-dir",
        default=os.environ.get("SCENIC_DIR", "results/SCENIC"),
        help="Directory containing HPAP_V3_regulons.csv",
    )
    args = parser.parse_args()

    scenic_dir = args.scenic_dir
    input_file = os.path.join(scenic_dir, "HPAP_V3_regulons.csv")
    output_file = os.path.join(scenic_dir, "HPAP_V3_regulon_targets.tsv")

    regulons = load_signatures(input_file)

    with open(output_file, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["regulon", "TF", "target", "weight"])

        for regulon in regulons:
            tf = getattr(
                regulon,
                "transcription_factor",
                regulon.name.split("(")[0],
            )
            for target, weight in regulon.gene2weight.items():
                writer.writerow([regulon.name, tf, target, weight])

    print("Number of regulons:", len(regulons))
    print("Saved:", output_file)


if __name__ == "__main__":
    main()
