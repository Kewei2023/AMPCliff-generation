# AMPCliff: quantitative definition and benchmarking of activity cliffs in antimicrobial peptides

Activity cliff (AC) is a phenomenon that two molecules differ by a small structural alternation but exhibit significantly changed biochemical activities. The AC in small molecules has been extensively investigated but limited knowledge is accumulated about the AC phenomenon in peptides with canonical amino acids. This study introduces a quantitative definition and benchmarking framework AMPCliff for the AC phenomenon in antimicrobial peptides (AMPs) composed by canonical amino acids. A comprehensive analysis of the existing AMP datasets reveals the significant prevalence of AC within AMPs. AMPCliff quantifies the AC between two peptides with the normalized BLOSUM62 similarity score in their alignment at least 0.9, and two-fold or larger minimum inhibitory concentration (MIC) change. This study establishes a benchmark dataset of paired AMPs in Staphylococcus aureus from the publicly available AMP dataset GRAMPA, and conducts a rigorous procedure to evaluate various AMP AC prediction models, including nine machine learning, four deep learning algorithms, four masked language models, and four generative language models. The experimental data reveals that these models are capable of detecting AMP AC events and the pre-trained protein language ESM2 model demonstrates superior performance across the comparative evaluations. The predictive model based on ESM2 with 33 layers only achieves the Spearman correlation coefficient (SCC) 0.50 for the regression task of the MIC values on the benchmark dataset. So the predictive performance of AMP AC remains to be further improved. Source code and additional resources are available at https://www.healthinformaticslab.org/supp/ or https://github.com/Kewei2023/AMPCliff-generation. The benchmarking code can be found at https://github.com/Kewei2023/AMPCliff.


## Versions

Kewei Li (kwbb1997@gmail.com) or Fengfeng Zhou (FengfengZhou@gmail.com)

Version 1.00


## Table of Contents

- [Installation](#installation)
- [Usage](#usage)

- [Contact](#contact)

## Installation

```
git clone git@github.com:Kewei2023/AMP-Cliff-generation.git

cd AMPCliff-generation

conda create -f environment.yaml
conda activate AMPCliff
```

## Usage
### Generate AMP-Cliff pairs
- **condition:** `blosum62 average` by default, supporting `levenstein aligned`,`levenstein`,`blosum62 average`,`tanimoto average`,`sequence identity`,`all`.


- **diff:** dilution difference,`5` by default, should be `int` like: `2`,`3`,`4`,`5`

- **threhsold:** similarity threshold for `blosum62 average` and `tanimoto average`,`0.9` by default, should be `float` like: `0.85`,`0.9`,`0.95`
```
python generate_cliffs.py -d xxx.csv -c "blosum62 average" -f 5 -t 0.9
```


### Data partition

```
python data_partition.py -d xxx.csv -p xxx.csv -t 0.9
```

the results will be saved in `./data/`

## Contact 
please contact kwbb1997@gmail.com for any question

## Update

2024-07-11

