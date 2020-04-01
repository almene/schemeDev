import vcfpy
from argparse import ArgumentParser
from Bio import SeqIO
import os
import csv


def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser_inner = ArgumentParser(
        description='BioHansel Scheme development accessory script')
    parser_inner.add_argument('--in_vcf', type=str, required=True, help='vcf File of SNPs')
    parser_inner.add_argument('--outdir', type=str, required=False, default=os.getcwd(),
                              help='Output Directory to put results')
    parser_inner.add_argument('--keep_csv', type=str, required=False,
                              help='csv file with comma deliminated positions to be kept for the reduced fasta file')
    parser_inner.add_argument('--remove_csv', type=str, required=False,
                              help='csv file with comma deliminated positions to be removed for the reduced fasta file')
    parser_inner.add_argument('--o', type=str, required=False, default="Reduced",
                              help='name for output files')
    return parser_inner.parse_args()


def main():
    args = parse_args()
    vcf_file = args.in_vcf
    out = args.o
    out_path = args.outdir
    try:
        reader = vcfpy.Reader.from_path(vcf_file)
    except FileNotFoundError:
        print(
            f"There was an error in reading {vcf_file}.  "
            f"Verify that the file exists")
        raise SystemExit(0)
    except vcfpy.exceptions.IncorrectVCFFormat:
        print(
            f"There was an error in reading {vcf_file}.  "
            f"Verify that the file is in the correct format")
        raise SystemExit(0)
    if not args.keep_csv:
        infile = args.remove_csv
        Opp = 0
    elif not args.remove_csv:
        infile = args.keep_csv
        Opp = 1
    else:
        print("Select one reduction method: either list positions to keep or list positions to remove")
        raise SystemExit(0)
    try:
        with open(infile) as tsvfile:
            keep = csv.reader(tsvfile, delimiter=",")
            for record in keep:
                pos_list = record
    except FileNotFoundError:
        print("There was an error reading the csv file.  Check the file name and try again")
        raise SystemExit(0)

    fasta_dict = dict()
    samples = reader.header.samples.names
    writer = vcfpy.Writer.from_path(f'{out}.vcf', reader.header)
    for record in reader:
        # skip over the metadata lines at the top of the file
        if not record.is_snv():
            continue
        position = f"{record.POS}"
        if Opp == 1 :
            if position in pos_list:
                writer.write_record(record)
                # pull the positional information and the reference value and all the alt values from the vcf file
                line = [record.CHROM, record.POS, record.REF]
                line += [alt.value for alt in record.ALT]

                # initialize variables and pull the ref and alt SNPs
                alt_base = record.ALT[0].value
                ref_base = record.REF
                tracker = 0

                # go through the samples and divide them by alt vs ref bases
                for call in record.calls:
                    state = int(call.data.get('GT'))
                    sample_name = samples[tracker]
                    if sample_name not in fasta_dict.keys():
                        fasta_dict[sample_name] = ""
                    if state == 0:
                        fasta_dict[sample_name] += ref_base
                    else:
                        fasta_dict[sample_name] += alt_base
                    tracker += 1
        if Opp == 0 :
            if position not in pos_list:
                writer.write_record(record)
                # pull the positional information and the reference value and all the alt values from the vcf file
                line = [record.CHROM, record.POS, record.REF]
                line += [alt.value for alt in record.ALT]

                # initialize variables and pull the ref and alt SNPs
                alt_base = record.ALT[0].value
                ref_base = record.REF
                tracker = 0

                # go through the samples and divide them by alt vs ref bases
                for call in record.calls:
                    state = int(call.data.get('GT'))
                    sample_name = samples[tracker]
                    if sample_name not in fasta_dict.keys():
                        fasta_dict[sample_name] = ""
                    if state == 0:
                        fasta_dict[sample_name] += ref_base
                    else:
                        fasta_dict[sample_name] += alt_base
                    tracker += 1

    fasta_file = open(os.path.join(out_path, f"{out}.fasta"), "w+")
    for key in fasta_dict.keys():
        fasta_file.write(f">{key}\n")
        fasta_file.write(f"{fasta_dict[key]}\n")


if __name__ == '__main__':
    main()
