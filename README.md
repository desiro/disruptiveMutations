# [<samp>disruptiveMutations</samp>](https://github.com/desiro/disruptiveMutations)
[![License: GPL v3](https://img.shields.io/badge/License-GPL_v3-bd0000.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python v3.9.7](https://img.shields.io/badge/Language-Python_v3-75a8d3.svg)](https://www.python.org/)
[![Conda v4.11.0](https://img.shields.io/badge/Uses-Conda-43b02a.svg)](https://docs.conda.io/en/latest/miniconda.html)

***

## Description

This tool reduces the interaction strength between a particular area of interest in the first sequence in relation to all other sequences. Ideally, the tool interrupts any potential interaction between this area and all other sequences. Individual predictions are performed with the RNAcofold python site-package of the ViennaRNA Package 2.5.0.

### Mandatory Prerequisites

* [![Python v3.9.7](https://img.shields.io/badge/Python_v3.9.7-75a8d3.svg)](https://www.python.org/downloads/release/python-397/)
* [![ViennaRNA v2.5.0](https://img.shields.io/badge/ViennaRNA_v2.5.0-006795.svg)](https://www.tbi.univie.ac.at/RNA/)

### Optional Prerequisites

* [![Conda v4.11.0](https://img.shields.io/badge/Conda_v4.11.0-43b02a.svg)](https://docs.conda.io/en/latest/miniconda.html)

***

## Installation

To run <samp>disruptiveMutations</samp>, I recommend using Miniconda and following the steps below. If this is the first time using conda, you should probably restart your shell after the installation of Miniconda. The following will demonstrate the installation and set up of Miniconda on Linux, which should be similar on other platforms. For Windows 10 users, I advise using the [Ubuntu 20.04 LTS](https://www.microsoft.com/en-us/p/ubuntu-2004-lts/9n6svws3rx71?cid=msft_web_chart) subsystem. More information can be found on the [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and [Bioconda](https://bioconda.github.io/user/install.html) pages.

### Conda Installation

Installing Miniconda:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Updating Miniconda and setting channels:
```
conda update conda
conda update python
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Installing Conda packages:
```
conda create --name disruptiveMutations python=3.9.7
conda activate disruptiveMutations
conda install -c bioconda viennarna=2.5.0
git clone https://github.com/desiro/disruptiveMutations.git
cd disruptiveMutations
```

### Installation without Conda

Installing the ViennaRNA package on Linux:
```
tar -zxvf ViennaRNA-2.5.0.tar.gz
cd ViennaRNA-2.5.0
./configure --with-python3
make
sudo make install
```

Installing the ViennaRNA package on MAC:
```
tar -zxvf ViennaRNA-2.5.0.tar.gz
cd ViennaRNA-2.5.0
./configure --enable-universal-binary --with-python3
make
sudo make install
```

***

## Examples

The basic input for ```disruptiveMutations.py``` includes the following parameters:
* ```-pfx``` - name of the output folder
* ```-fsa``` - input fasta file

### Fasta Formation

The tool always takes the first sequence in the fasta file as the target mutation sequence. If there is no specified range, the whole sequence will be taken and compared against all other sequences present in the fasta file. 

### Basic Example

```
python disruptiveMutations.py -pfx example -fsa example.fa -pss 32 -pse 96 -thr 4 -ovr
```

### Options

For more command line options, see the [manual](https://github.com/desiro/disruptiveMutations/blob/master/manual.md).

***

## Authors

* [Daniel Desirò](https://github.com/desiro)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Reference

Please cite <samp>DisruptiveMutations</samp> if you find our tool useful.

```
D. Desirò, A. Borodavka, and M. Marz.
"DisruptiveMutations: disrupting functional long-range RNA-RNA interactions in RNA viruses."
In Preparation, 2022.
```
