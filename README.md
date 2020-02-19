# bioCanon

Biohansel requires a strain specific fasta format scheme to be able to classify samples.
Schemes exist for some common clonal pathogens but there are many untapped applications.
This module includes functions to aid in the generation of k-mers linked to predefined groups
to serve as the basis for automatic scheme development

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for 
development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites
Several python packages are required to run bioCanon:
vcfpy
argparse
biopython
ete3
numpy
pandas

All of these packages can be installed via pip
```
pip install vcfpy
```
Note: vcfpy has several python dependencies that will also need to be installed (cython, pysam) some of these 
require some libraries to be installed via apt-get (zlib1g-dev,  libbz2-dev and liblzma-dev on Ubuntu)

### Installing

After all dependencies have been installed bioCanon can be installed via pip from the github repository

```
pip install git+https://github.com/almene/schemeDev.git#egg=schemeDev
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Tests for the module can be found in the respository under the tests subfolder.  They are pytest compatible scripts and once downloaded can be read using the python -m pytest command

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Amanda Saunders** - *Initial work* - [Almene](https://github.com/almene)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* James Robertson and Justin Schonfeld for their help in development
