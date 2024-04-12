# AMP-Cliff-generation

A brief description of what this project does and who it's for.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)

- [Contact](#contact)

## Installation

```
git clone git@github.com:Kewei2023/AMP-Cliff-generation.git

cd AMP-Cliff-generation

conda create -f environment.yaml
conda activate AMPCliff
```

## Usage
### Generate AMP-Cliff pairs
**condition:** `blosum62 average` by default, supporting `levenstein aligned`,`levenstein`,`blosum62 average`,`tanimoto average`,`sequence identity`,`all`.


**diff:** dilution difference,`5` by default, should be int like: `2`,`3`,`4`,`5`
```
python generate_cliffs.py -d xxx.csv -c "blosum62 average" -f 5
```


### Data partition

```
python data_partition.py -d xxx.csv -p xxx.csv
```

the results will be saved in `./data/`

## Contact 
please contact kwbb1997@gmail.com for any question