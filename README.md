[<img src="img/rice_pine.png" width="1000"/>]()

# genomic prediction for GBLUP estimated breeding values (EBVs) associated to rice and pine traits

### 🎯 Objective

This repository contains R scripts designed for reproducible data analysis and results, aligned with the FAIR principles. The scripts perform data reformatting and GBLUP EBVs computation to be used as phenotypes for genomic prediction.

### 💻 Instructions

Download the ```double_gblup_analyses``` repository in the current user's directory on a computing cluster or personal computer using one of the following commands :

  *  ```git clone git@github.com:ljacquin/double_gblup_analyses.git``` <p> </p>
    or
  * ```git clone https://github.com/ljacquin/double_gblup_analyses.git``` 
  <p> </p>
  
  ⚠️ Make sure``` git``` is installed beforehand; if not, install it with ```sudo apt install git```.
  <p> </p>

* Given that ```R ≥ 4.1.2``` is already installed, within the ```double_gblup_analyses``` folder use the following command to install and test ```double_gblup_analyses``` required ```R``` libraries : 

  * ```R -q --vanilla < src/requirements.R```
  * ```R -q --vanilla < src/test_requirements.R```
  <p> </p>
  
* The ```R``` scripts ```*data_reformatting.R```, ```*spat_hetero_correct_per_env_trait.R``` and ```*gblup_ebv_computation.R``` in the ```double_gblup_*/src/data_treatment_and_analysis/``` folders 
perform 0) data reformatting, 1) spatial heterogeneity correction, and 2) GBLUP EBVs computation.

* The ```R``` script ```double_gblup_*/src/genomic_prediction_and_analysis/genomic_prediction.R``` performs, for each trait, the genomic prediction tasks and analyses for the GBLUP EBVs

⚠️ The tasks and analyses performed by the ```R``` scripts in the ```double_gblup_analyses``` repository can be run in either ```Unix/Linux``` or ```Windows``` environments, as long as ```R``` and the necessary libraries are installed. For local computations in ```RStudio```, ensure that the ```computation_mode``` variable is set to "local" in the ```R``` scripts located in ```src/```.


