#!/usr/bin/env python
# Changes: python3, argparse, gff3
# check if python >= 3.6

import argparse
import sys


def main(input_file, output_file):
    """Converts a file from the BED (Browser Extensible Data) format to the GFF3 (General Feature Format version 3) format."""
    skipped_lines = 0
    first_skipped_line = 0
    out = open(output_file, 'w')
    out.write('##gff-version 3\n')
    out.write('##bed_to_gff3_converter.py\n\n')
    i = 0
    prev_identifier = ""
    for i, line in enumerate(open(input_file)):
        complete_bed = False
        line = line.rstrip('\r\n')
        if line and not line.startswith('#') and not line.startswith('track') and not line.startswith('browser'):
            try:
                elems = line.split('\t')
                if len(elems) == 12:
                    complete_bed = True
                chrom = elems[0]
                if complete_bed:
                    feature = 'mRNA'
                else:
                    try:
                        feature = elems[3]
                    except Exception:
                        feature = f'feature{i + 1}'
                start = int(elems[1]) + 1
                end = int(elems[2])
                try:
                    score = elems[4]
                except Exception:
                    score = '0'
                try:
                    strand = elems[5]
                except Exception:
                    strand = '+'
                try:
                    identifier = elems[3]
                    if identifier == prev_identifier.split(".")[0]:
                        path_number = int(prev_identifier.split(".")[1].lstrip("mrna")) + 1
                        identifier = f"{identifier}.mrna{path_number}"
                    else:
                        identifier = f"{identifier}.mrna1"
                    prev_identifier = identifier
                except Exception:
                    identifier = f'ID.{i + 1}'
                out.write(f'{chrom}\tbed_to_gff3_converter\t{feature}\t{start}\t{end}\t{score}\t{strand}\t.\tID={identifier}\n')
                if complete_bed:
                    block_count = int(elems[9])
                    block_sizes = elems[10].split(',')
                    block_starts = elems[11].split(',')
                    for j in range(block_count):
                        exon_start = int(start) + int(block_starts[j])
                        exon_end = exon_start + int(block_sizes[j]) - 1
                        out.write(f'{chrom}\tbed_to_gff3_converter\texon\t{exon_start}\t{exon_end}\t{score}\t{strand}\t.\tParent={identifier}\n')
            except Exception:
                skipped_lines += 1
                if not first_skipped_line:
                    first_skipped_line = i + 1
        else:
            skipped_lines += 1
            if not first_skipped_line:
                first_skipped_line = i + 1
    out.close()
    info_msg = f'{i + 1 - skipped_lines} lines converted to GFF3.'
    if skipped_lines > 0:
        info_msg += f'\nSkipped {skipped_lines} blank/comment/invalid lines starting with line #{first_skipped_line}.'
    print(info_msg)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converts a file from the BED (Browser Extensible Data) format to the GFF3 (General Feature Format version 3) format.')
    parser.add_argument('input_file', help='Input file in the BED format.')
    parser.add_argument('output_file', help='Output file in the GFF3 format.')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    main(**vars(args))