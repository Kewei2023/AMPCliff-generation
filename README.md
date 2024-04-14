# AMP-Cliff-generation

Activity cliff (AC) is a phenomenon that pairs of similar compounds that only differ by a small structural modification but exhibit a large difference in their binding affinity for a given target. The AC of small molecules has been extensively investigated but limited knowledge is accumulated about the AC phenomenon in peptides with canonical amino acids. This study introduces a quantitative definition and benchmarking framework AMPCliff for the AC phenomenon in antimicrobial peptides (AMPs). A comprehensive analysis of the existing AMP dataset reveals a significant prevalence of AC within AMPs. AMPCliff quantifies the activities of AMPs by the metric minimum inhibitory concentration (MIC), and defines 0.9 as the minimum threshold for the normalized BLOSUM62 similarity score between a pair of aligned peptides with at least two-fold MIC changes. This study establishes a benchmark dataset of paired AMPs in Staphylococcus aureus from the publicly available AMP dataset GRAMPA, and conducts a rigorous procedure to evaluate various AMP AC prediction models, including nine machine learning, four deep learning algorithms, four masked language models, and four generative language models. Our analysis reveals that these models are capable of detecting AMP AC events and the pre-trained protein language ESM2 model demonstrates superior performance across the evaluations. The predictive performance of AMP activity cliffs remains to be further improved, considering that ESM2 only achieves the Spearman=0.50 of AC with 2-fold change. Source code and additional resources are available at https://github.com/Kewei2023/AMPCliff-generation

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