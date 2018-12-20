This repository consists of jupyter notebooks with python, bash and R code for reproducing results of paper Zenin et al. "Identification of 12 genetic loci associated with human healthspan".

Requirements:

- Any suitable linux OS, (e.g. Ubuntu 16.04 LTS)
- ascp version 3.7.4.147133 https://downloads.asperasoft.com/en/downloads/50
- aws-cli/1.11.13 Python/3.5.2 Linux/4.4.0-1070-aws botocore/1.4.70 https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-install.html
- Ega Demo Download Client  Version: 2.2.2 https://ega-archive.org/download/using-ega-download-client
- bgenix version: 1.2-dev https://bitbucket.org/gavinband/bgen/wiki/bgenix
- cat-bgen https://bitbucket.org/gavinband/bgen/wiki/cat-bgen
- PLINK v1.90b6.5 64-bit (13 Sep 2018) https://www.cog-genomics.org/plink/1.9/
- PLINK v2.00a2LM AVX2 Intel (22 Sep 2018) https://www.cog-genomics.org/plink/2.0/
- PAINTOR v3.0 https://github.com/gkichaev/PAINTOR_V3.0
- Genome-wide Complex Trait Analysis (GCTA) * version 1.91.1 beta http://cnsgenomics.com/software/gcta/

All python dependencies are listed at requirements.txt file.

Please note, this work was done using very large datasets from UK Biobank, therefore it is not feasable to run this analysis on regular desktop or laptor computers, some high-performance hardware servers needed, such as provoded by Amazon Web Services, Google Cloud Platform or any other IaaS providers available at pay-per-use basis.

To store and process data you will need at least 5 TB hard drive space, 60-120 GB RAM, and 8-64x cores high-performance CPU. The full run starting from downloading data and further pre-processing and running GWAS may take up to two weeks to complete.

Content:

01_Genetics_data_processing.ipynb - Downloading UK Biobank data and preparing reference datasets
02_1_Cox-Gomperz_healthspan_model.ipynb - Defining healthspan and calculating it from UK Biobank data, fitting Cox-Goperz model for covariates, preparing datasets for replication.
03_Functional_annotation.ipynb - Genetic functional annotation pipelines for PAINTOR, COJO and VEP
04_Genetic_correlations.ipynb - Computing genetics correlations for calls data.

Instructions:

The main results can be obtained by running gwas_on_calls.py or gwas_on_imputed.py python scripts, all demo data required for these runs available at data directory.