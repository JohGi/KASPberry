#!/usr/bin/env python3

import argparse
from ids import make_snp_id

FAILURE_REASON = "non_diagnostic_allele_pattern"


def read_list(path):
    """Read one value per line."""
    with open(path, encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def parse_header(line):
    """Return column index mapping."""
    header = line.rstrip("\n").split("\t")
    return {name: i for i, name in enumerate(header)}


def is_discriminant(fields, idx, group_a, group_b):
    """Return True if SNP discriminates the two groups."""
    vals_a = [fields[idx[sample]] for sample in group_a]
    vals_b = [fields[idx[sample]] for sample in group_b]

    return (
        len(set(vals_a)) == 1
        and len(set(vals_b)) == 1
        and vals_a[0] != vals_b[0]
    )


def write_status_row(fout, fields, passed):
    """Write one SNP filtering status row."""
    block_id = fields[0].removesuffix(".aln")
    aln_pos = fields[1]
    snp_id = make_snp_id(block_id, aln_pos)

    status = "PASS" if passed else "FAIL"
    failure_reason = "" if passed else FAILURE_REASON

    fout.write(
        f"{snp_id}\t{block_id}\t{aln_pos}\t"
        f"{status}\t{failure_reason}\n"
    )


def write_snp_filter_status(input_path, output_path, group_a, group_b):
    """Evaluate all SNPs and write their diagnostic filtering status."""
    with (
        open(input_path, encoding="utf-8") as fin,
        open(output_path, "w", encoding="utf-8") as fout,
    ):
        fout.write(
            "snp_id\tblock_id\taln_pos\tstatus\tfailure_reason\n"
        )

        idx = None

        for line in fin:
            if line.startswith("#CHROM"):
                idx = parse_header(line)
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            passed = (
                is_discriminant(fields, idx, group_a, group_b)
                if group_a
                else True
            )

            write_status_row(fout, fields, passed)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--group-a-file", required=True)
    parser.add_argument("--group-b-file", required=True)
    args = parser.parse_args()

    group_a = read_list(args.group_a_file)
    group_b = read_list(args.group_b_file)

    write_snp_filter_status(
        args.input,
        args.output,
        group_a,
        group_b,
    )


if __name__ == "__main__":
    main()
