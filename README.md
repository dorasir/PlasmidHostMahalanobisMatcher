# PlasmidHostMahalanobisMatcher

PlasmidHostMahalanobisMatcher is a software for predicting the possible host for given plasmids based on multiple features.

## Prerequisites

### Program
* This software requires Python 3, install [Python3](https://www.python.org/downloads/release/python-363/) or [Anaconda3](https://www.anaconda.com/download/)

### Required packages
* The packages used in this software is included in ```requirements.txt```
* If you are using ```pip```, run ```pip install -r requirements.txt``` in the console
* If you are using ```conda```, run ```conda install --file requirements.txt``` in the console


### Installation
* Download the repository and unzip or clone the repository 
  ```
  git clone https://github.com/dorasir/PlasmidHostMahalanobisMatcher
  ```
* For MacOS, run the following command to setup
  ```
  MACOSX_DEPLOYMENT_TARGET=10.9 CC=g++ python setup.py install --install-platlib=src
  ```
* For Linux, run the following command to setup
  ```
  CC=g++ python setup.py install --install-platlib=src
  ```

## Data availability
Our constructed host bacteria database is available for this software, the database can be downloaded at [GoogleDrive](https://drive.google.com/file/d/13VkJF5Bcw3Y2CYhFMOuesil1eVKvndCZ/view?usp=sharing). The content of the downloaded archive should then be extracted to ```data``` folder of the software directory.

## Usage

```
python PlasmidHostMahalanobisMatcher.py [-h] -q QUERY_PLASMID_DIR -o OUTPUT [-t NUM_THREADS] [-n topN] [-i INTERMEDIATE_DIR]
```

### Options
```
-h, --help            show this help message and exit
-q QUERY_PLASMID_DIR  Directory containing query plasmid genomes with .fasta or .fa suffix
-o OUTPUT             Output file
-t NUM_THREADS        Number of threads to use. Default = 1
-n topN               Number of top predictions written to the output files. All predictions will be output if there is a tie in score. Default = 1
-i INTERMEDIATE_DIR   Directory storing intermediate result. Default = ./intermediate_res
```

### Examples

To predict the host of plasmids of interest, put the genome of the plasmids into a dedicated folder, for example, ```plasmids```. The prediction can then be made by running 
```
python PlasmidHostMahalanobisMatcher.py -q plasmids
```
By default the prediction result will be stored in ```output.txt``` under the software folder.