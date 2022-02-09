# disruptiveMutations

```disruptiveMutations``` reduces the interaction strength between a given area of the first sequence in regard to all other sequences. The tool is written in ```Python 3.7.1``` and the calculations are performed with the ```RNAcofold``` python site-package of the ```ViennaRNA Package 2.4.13```. For command-line options, please refer to the [manual](https://github.com/desiro/disruptiveMutations/blob/master/manual.md)

## Mandatory Prerequisites

* [python 3.7.1](https://www.python.org/downloads/release/python-385/)
* [viennaRNA 2.4.13](https://www.tbi.univie.ac.at/RNA/documentation.html#install)

## Optional Prerequisites

* [Miniconda3](https://docs.conda.io/en/latest/miniconda)

## Installation with Miniconda

Please download the [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) application for your system. The following will demonstrate the installation and set up of miniconda on Linux, which should be similar on other platforms.

```
bash Miniconda3-latest-Linux-x86_64.sh -p ~/miniconda3
conda create --name disruptiveMutations
conda activate disruptiveMutations
conda install -c bioconda viennarna
```

## Examples

The basic input for ```disruptiveMutations.py``` includes the following parameters:
* ```-pfx``` - name of the output folder
* ```-fsa``` - input fasta file

### Fasta Formation

The tool always takes the first sequence in the fasta file as the target mutation sequence. If there is no specified range, the whole sequence will be taken and compared against all other sequences present in the fasta file. 

### Basic Example

```
python3 disruptiveMutations.py -pfx example -fsa example.fa -pss 32 -pse 96 -thr 4 -ovr 
```

## Authors

* [Daniel Desirò](https://github.com/desiro)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Reference

Please cite [rotavirusPaper](https://doi.org/10.1101/424002) if you find our tool useful.

