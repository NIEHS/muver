# muver
muver is designed to call mutations in outgrowths from mutation accummulation experiments. Developed to address specific experimental challenges in yeast, muver is designed to leverage high sequencing depth and access to data from both ancestor and progeny.

## Requirements
muver was developed using Python 2.7.13. In addition to requirements specified in setup.py, muver requires installation of the following tools:
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [GATK, version 3.7-0](https://software.broadinstitute.org/gatk/download/)
* [picard](https://broadinstitute.github.io/picard/)
* [samtools](http://www.htslib.org/download/)

## Installation
Proper function of muver requires the paths to depencies to be set.  To do this, manually set the paths in `paths.cfg` using a text editor.

After the correct paths have been set, install muver with the following command:
```
python setup.py install
```

## Usage
All of muvers functions may be accessed using its command line interface. General usage is as follows:
```
muver COMMAND [OPTIONS] [ARGS]...
```
A list of commands can be found by using the following:
```
muver --help
```
Details about each command can be found by using the following:
```
muver COMMAND --help
```
See the [manual](docs/manual.md) for further usage details.

## Authors
muver was conceptualized by Scott Lujan and Adam Burkholder. muver was written by Adam Burkholder and Christopher Lavender.

## License
This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgments
This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [audreyr/cookiecutter-pypackage](https://github.com/audreyr/cookiecutter-pypackage) project template.
