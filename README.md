# bioCanon

Biohansel requires a strain-specific fasta format scheme to be able to classify samples.
Schemes exist for some common clonal pathogens, but there are many untapped applications.
This module includes functions to aid in the generation of k-mers linked to predefined groups
to serve as the basis for automatic scheme development.  The only test of k-mer quality is the
exclusion of k-mers with degenerate bases.  Further filtering of output kmers will be required
to improve the quality of the biohansel scheme.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for 
development and testing purposes.

### Prerequisites
Several python packages are required to run bioCanon:
vcfpy,
argparse,
biopython,
ete3,
numpy, and
pandas

All of these packages can be installed via pip
```
pip install vcfpy
pip install argparse
pip install biopython
pip install ete3
pip install numpy
pip install pandas
```
Note: vcfpy has several python dependencies that will also need to be installed (cython, pysam) some of these 
require some libraries to be installed via apt-get (zlib1g-dev,  libbz2-dev and liblzma-dev on Ubuntu)

### Installing

After all dependencies have been installed bioCanon can be installed via pip from the GitHub repository

```
pip install git+https://github.com/almene/schemeDev.git#egg=biocanon-almene
```

### Using the Module

Once it has been installed, the module can be run from the command line.  There are several arguments that the module accepts.  Some of them are required and some are optional.

```
python3 -m bioCanon --in_vcf tests/examples/testing.vcf --in_nwk tests/examples/testing.nwk --reference tests/examples/ref.fasta --min_snps 2 --min_members 5 --min_parent 5 --outdir testcases
```

#### Arguments

##### --in_vcf
* Required

vcfpy reader compatible vcf file.  See vcfpy documentation for specific requirements

##### --in_nwk or --group_info
* Either --in_nwk or --group_info required

--in_nwk requires a Newick format tree file
```
((((((((((((Q,R),(O,P)),J),I),H),G),(S,F)),B),(C,(T,U))),(((M,N),(K,L)),E)),(D,V)),A);
```
--group_info requires a tsv format file that contains the group information

```
A    1    
B    2    2.1    2.1.1    2.1.1.1    2.1.1.1.1
C    2    2.1    2.1.1    2.1.1.2    2.1.1.2.1
D    2    2.2    
E    2    2.1    2.1.2    2.1.2.1    2.1.1.1.2
F    2    2.1    2.1.1    2.1.1.1    2.1.1.1.2    2.1.1.1.2.2
G    2    2.1    2.1.1    2.1.1.1    2.1.1.1.2    2.1.1.1.2.1    2.1.1.1.2.1.2
H    2    2.1    2.1.1    2.1.1.1    2.1.1.1.2    2.1.1.1.2.1    2.1.1.1.2.1.1    2.1.1.1.2.1.1.2
I    2    2.1    2.1.1    2.1.1.1    2.1.1.1.2    2.1.1.1.2.1    2.1.1.1.2.1.1    2.1.1.1.2.1.1.1    2.1.1.1.2.1.1.1.2
J    2    2.1    2.1.1    2.1.1.1    2.1.1.1.2    2.1.1.1.2.1    2.1.1.1.2.1.1    2.1.1.1.2.1.1.1    2.1.1.1.2.1.1.1.1    2.1.1.1.2.1.1.1.1.2
K    2    2.1    2.1.2    2.1.2.2    2.1.2.2.2
L    2    2.1    2.1.2    2.1.2.2    2.1.2.2.2
M    2    2.1    2.1.2    2.1.2.2    2.1.2.2.1
N    2    2.1    2.1.2    2.1.2.2    2.1.2.2.1
O    2    2.1    2.1.1    2.1.1.1    2.1.1.1.2    2.1.1.1.2.1    2.1.1.1.2.1.1    2.1.1.1.2.1.1.1    2.1.1.1.2.1.1.1.1    2.1.1.1.2.1.1.1.1.1    2.1.1.1.2.1.1.1.1.1.2
P    2    2.1    2.1.1    2.1.1.1    2.1.1.1.2    2.1.1.1.2.1    2.1.1.1.2.1.1    2.1.1.1.2.1.1.1    2.1.1.1.2.1.1.1.1    2.1.1.1.2.1.1.1.1.1    2.1.1.1.2.1.1.1.1.1.2
Q    2    2.1    2.1.1    2.1.1.1    2.1.1.1.2    2.1.1.1.2.1    2.1.1.1.2.1.1    2.1.1.1.2.1.1.1    2.1.1.1.2.1.1.1.1    2.1.1.1.2.1.1.1.1.1    2.1.1.1.2.1.1.1.1.1.1
R    2    2.1    2.1.1    2.1.1.1    2.1.1.1.2    2.1.1.1.2.1    2.1.1.1.2.1.1    2.1.1.1.2.1.1.1    2.1.1.1.2.1.1.1.1    2.1.1.1.2.1.1.1.1.1    2.1.1.1.2.1.1.1.1.1.1
S    2    2.1    2.1.1    2.1.1.1    2.1.1.1.2    2.1.1.1.2.2
T    2    2.1    2.1.1    2.1.1.2    2.1.1.2.2
U    2    2.1    2.1.1    2.1.1.2    2.1.1.2.2
V    2    2.2
```
Above is an example of the required format

##### --reference
* Required

fasta format reference genome file

##### --min_snps, --min_members, and --min_parent
* Optional

Integer variables that restrict what is considered a valid group: minimum number of snps that support a division, minimum members required to generate a valid group and minimum size difference between a subgroup and its parent group required for the subgroup to be considered valid

* Defaults : 2, 5, 2
##### --flanking
* Optional

Integer number of bases flanking the SNP site determines the k-mer size

* Default: 15

##### --outdir
* Optional

Name for an output directory

### Outputs
Three files will be produced by a successful run of the bioCanon module:
* codes.log
  * This file contains the samples used to build the scheme along with the biohansel code that they were assigned by the scheme
  ```
  Q	1.6
  R	1.6
  O	1.6
  P	1.6
  J	1.6
  I	1
  H	1
  G	1.5
  S	1.2
  ...
  ```
* a S[x]G[y]_biohansel.fasta file
  * This file contains all the non-degenerate positve and negative k-mers required for a biohansel scheme.  This will require further filtering for optimal scheme development
  * Naming convention is to include parameter information [x] is the required snp support per valid group and [y] is the minimum group size 
    ```
    >305-1
    GAACAGTGCGGGCTTTTTTTTCGACCAGAGA
    >negative305-1
    GAACAGTGCGGGCTTGTTTTTCGACCAGAGA
    >2685-1
    GGCGTTCTACAGCCATTATTATCAGCCCTTG
    >negative2685-1
    GGCGTTCTACAGCCACTATTATCAGCCCTTG
    >5021-1
    ATGACCCGCCAGTAATCATTGCGCCCGGTGG
    ...
    ```
* a S[x]G[y]_biohansel.log file
  * This file contains a more human-readable version of the information contained in the biohansel.fasta file
  * Naming convention is to include parameter information [x] is the required snp support per valid group and [y] is the minimum group size
    ```
    305	1	GAACAGTGCGGGCTTTTTTTTCGACCAGAGA
    negative305	1	GAACAGTGCGGGCTTGTTTTTCGACCAGAGA
    2685	1	GGCGTTCTACAGCCATTATTATCAGCCCTTG
    negative2685	1	GGCGTTCTACAGCCACTATTATCAGCCCTTG
    5021	1	ATGACCCGCCAGTAATCATTGCGCCCGGTGG
    negative5021	1	ATGACCCGCCAGTAACCATTGCGCCCGGTGG
    5392	1	CAGATTCACCACCACACGATCGCCCTGCGCC
    negative5392	1	CAGATTCACCACCACCCGATCGCCCTGCGCC
    1180	1	CAGATCCCCTGTCTGTTTAAAAATACCGGTA
    negative1180	1	CAGATCCCCTGTCTGATTAAAAATACCGGTA
    ...
    ```
* a S[x]G[y]snp_report.txt file
  * This file contains a detailed report of the results of bioCanon for each snp in the vcf file
  * Naming convention is to include parameter information [x] is the required snp support per valid group and [y] is the minimum group size
    ```
    Position	Included	Supports	 Exclusion Reason
    1366	N		 no variability in the supplied vcf file
    68	N		 could not be matched to a branch point in the data.
    1017	N		 could not be matched to a branch point in the data.
    366	N		 a non-degenerate 31bp kmer could not be generated.
    945	N		 a non-degenerate 31bp kmer could not be generated.
    2050	N		 a non-degenerate 31bp kmer could not be generated.
    T6736	N		the number of size difference between its group and the parent group was less than 1.
    T6991	N		the number of size difference between its group and the parent group was less than 1.
    5110	Y	1.1
    9797	Y	1.6
    ...
    ```

## Running the tests

Tests for the module can be found in the repository under the tests subfolder.  They are pytest compatible scripts and once downloaded can be read using the python -m pytest command from the downloaded repository directory.

The current testing suite tests functions in main by supplying them with test data and comparing it against an expected result.

## Authors

* **Amanda Saunders** - *Initial work* - [Almene](https://github.com/almene)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* James Robertson and Justin Schonfeld for their help in the development
