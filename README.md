# AMP-Cliff-generation

This paper firstly defined the activity cliff phenomenon in antimicrobial peptides (AMP-Cliff), where structurally similar compounds exhibit significant variations in Minimum Inhibitory Concentrations (MIC) against identical bacterial strains. We compared the widely used conception of MMP-Cliff with AMP-Cliff, revealed a gap in the existing MMP-Cliff definition, which overlooks the structural similarity and evolutionary information of substructures within molecule pairs. Noting the synthetic nature of many molecules necessitating extended data collection, we capitalized on the prevalence of canonical peptides in published datasets to address this gap. Furthermore, we proposed a strict model evaluation framework to construct a thorough evaluation of recent algorithms for predicting AMP-Cliff, utilizing 9 machine learning and 4 deep learning approaches, 4 masked and 4 generative language models on peptides against S.aureus from the public AMP dataset GRAMPA. Our findings indicate that these models can discern low-frequency patterns to predict high-frequency signals, and ESM2 performed best in general. The code is available at https://github.com/Kewei2023/AMP-Cliff-generation

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