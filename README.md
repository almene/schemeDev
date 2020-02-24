# bioCanon

Biohansel requires a strain specific fasta format scheme to be able to classify samples.
Schemes exist for some common clonal pathogens but there are many untapped applications.
This module includes functions to aid in the generation of k-mers linked to predefined groups
to serve as the basis for automatic scheme development

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

After all dependencies have been installed bioCanon can be installed via pip from the github repository

```
pip install git+https://github.com/almene/schemeDev.git#egg=schemeDev
```

### Using the Module

Once it has been installed the module can be run from the comandline.  There are several arguments that the module accepts.  Some of them are required and some are optional.

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
A	1	
B	2	2.1	2.1.1	2.1.1.1	2.1.1.1.1
C	2	2.1	2.1.1	2.1.1.2	2.1.1.2.1
D	2	2.2	
E	2	2.1	2.1.2	2.1.2.1	2.1.1.1.2
F	2	2.1	2.1.1	2.1.1.1	2.1.1.1.2	2.1.1.1.2.2
G	2	2.1	2.1.1	2.1.1.1	2.1.1.1.2	2.1.1.1.2.1	2.1.1.1.2.1.2
H	2	2.1	2.1.1	2.1.1.1	2.1.1.1.2	2.1.1.1.2.1	2.1.1.1.2.1.1	2.1.1.1.2.1.1.2
I	2	2.1	2.1.1	2.1.1.1	2.1.1.1.2	2.1.1.1.2.1	2.1.1.1.2.1.1	2.1.1.1.2.1.1.1	2.1.1.1.2.1.1.1.2
J	2	2.1	2.1.1	2.1.1.1	2.1.1.1.2	2.1.1.1.2.1	2.1.1.1.2.1.1	2.1.1.1.2.1.1.1	2.1.1.1.2.1.1.1.1	2.1.1.1.2.1.1.1.1.2
K	2	2.1	2.1.2	2.1.2.2	2.1.2.2.2
L	2	2.1	2.1.2	2.1.2.2	2.1.2.2.2
M	2	2.1	2.1.2	2.1.2.2	2.1.2.2.1
N	2	2.1	2.1.2	2.1.2.2	2.1.2.2.1
O	2	2.1	2.1.1	2.1.1.1	2.1.1.1.2	2.1.1.1.2.1	2.1.1.1.2.1.1	2.1.1.1.2.1.1.1	2.1.1.1.2.1.1.1.1	2.1.1.1.2.1.1.1.1.1	2.1.1.1.2.1.1.1.1.1.2
P	2	2.1	2.1.1	2.1.1.1	2.1.1.1.2	2.1.1.1.2.1	2.1.1.1.2.1.1	2.1.1.1.2.1.1.1	2.1.1.1.2.1.1.1.1	2.1.1.1.2.1.1.1.1.1	2.1.1.1.2.1.1.1.1.1.2
Q	2	2.1	2.1.1	2.1.1.1	2.1.1.1.2	2.1.1.1.2.1	2.1.1.1.2.1.1	2.1.1.1.2.1.1.1	2.1.1.1.2.1.1.1.1	2.1.1.1.2.1.1.1.1.1	2.1.1.1.2.1.1.1.1.1.1
R	2	2.1	2.1.1	2.1.1.1	2.1.1.1.2	2.1.1.1.2.1	2.1.1.1.2.1.1	2.1.1.1.2.1.1.1	2.1.1.1.2.1.1.1.1	2.1.1.1.2.1.1.1.1.1	2.1.1.1.2.1.1.1.1.1.1
S	2	2.1	2.1.1	2.1.1.1	2.1.1.1.2	2.1.1.1.2.2
T	2	2.1	2.1.1	2.1.1.2	2.1.1.2.2
U	2	2.1	2.1.1	2.1.1.2	2.1.1.2.2
V	2	2.2
```
Above is an example of the required format

##### --reference
* Required

fasta format reference genome file

##### --min_snps, --min_members, and --min_parent
* Optional

Integer variables that restrict what is considered a valid group: minimum number of snps that support a division, miminum members required to generate a valid group and minimum size difference between a subgroup and its parent group required for the subgroup to be considered valid

* Defaults : 2, 5, 2
##### --flanking
* Optional

Integer number of bases flanking the SNP site, determines k-mer size

* Default: 15

##### --outdir
* Optional

Name for an output directory

### Outputs
Three files will be produced by a successful run of the bioCanon module:
* codes.log
  * contains the samples used to build the scheme along with the biohansel code that they were assigned by the scheme
* a biohansel.fasta file
  * a biohansel compatiple scheme file
* a buigabsek.log file
  * contains a more human readable version of the information contained in the biohansel.fasta file

Names for the biohansel.fasta and biohansel.log files are derived from the options used to construct them.  The naming structure is as follows:
```
S(min snps)G(mingroupsize)\_biohansel.fasta(or log)
```
## Running the tests

Tests for the module can be found in the respository under the tests subfolder.  They are pytest compatible scripts and once downloaded can be read using the python -m pytest command from the downloaded repository directory.

The current testing suite tests functions in main by supplying them with test data and comparing it against an expected result.

## Authors

* **Amanda Saunders** - *Initial work* - [Almene](https://github.com/almene)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* James Robertson and Justin Schonfeld for their help in development
