import vcfpy
from argparse import ArgumentParser
from Bio import SeqIO

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(
        description='BIO_Hansel Scheme development')
    parser.add_argument('--in_vcf', type=str, required=True, help='VCF File of SNPs')

    return parser.parse_args()

def main():
    args = parse_args()
    vcf_file = args.in_vcf
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
    fasta_dict = dict()
    samples = reader.header.samples.names
    for record in reader:
        # skip over the metadata lines at the top of the file
        if not record.is_snv():
            continue

        # pull the positional information and the reference value and all the alt values from the vcf file
        line = [record.CHROM, record.POS, record.REF]
        line += [alt.value for alt in record.ALT]

        # skip any vcf calls that have multiple alternative alleles
        if len(record.ALT) > 1:
            continue

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

    for key in fasta_dict.keys():
        print(f">{key}")
        print(fasta_dict[key])


if __name__ == '__main__':
    main()