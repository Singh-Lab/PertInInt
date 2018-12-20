# PertInInt

<div align="center"><img src="http://compbio.cs.princeton.edu/pertinint/PertInInt_workflow.png" alt="pipeline figure" title="PertInInt workflow" width="80%" /></div><br />

The goal of our analytical PertInInt method is to rapidly uncover proteins with significant enrichments of somatic **Pert**urbations **In** **Int**eraction and other functional sites. If you use data or scripts from this repository, please cite:

> Kobren, S.N., Chazelle, B. and Singh, M. (2019). "An integrative approach to identify preferentially altered interactions in human cancers." *Manuscript in preparation.*

### 1: Downloading preliminary data

* **Precomputed Tracks.** PertInInt uses *tracks* to model different functional regions of a protein. These tracks measure interaction interfaces, functional protein domains and conserved protein positions. To download the set of precomputed tracks used by PertInInt (which will be unzipped into a directory called `track_weights/`), run the following: 

  ```bash
  PERTININT_TRACKS="PertInInt-tracks_v0.tar.gz"
  wget http://compbio.cs.princeton.edu/pertinint/${PERTININT_TRACKS}
  tar -xvzf ${PERTININT_TRACKS}
  ```

* **Somatic Mutations.** PertInInt runs on any input .maf file. We tested PertInInt using a .maf file containing somatic mutations from all 33 TCGA cancer types, obtained from [NCI's Genomic Data Commons](https://gdc.cancer.gov) on December 6, 2018. To download this aggregated .maf file to test on, run the following: 

  ```bash
  if [ ! -d mafs ]; then mkdir mafs; fi
  AGGREGATE_CANCER="TCGA.Aggregate.muse.aggregated.somatic.maf.gz"
  wget http://compbio.cs.princeton.edu/pertinint/${AGGREGATE_CANCER} -O mafs/${AGGREGATE_CANCER}
  gzip -d mafs/${AGGREGATE_CANCER}
  ```

* **Expression Data.** We find that PertInInt's performance improves when limiting to those somatic mutations that fall into proteins that are expressed at TPM (transcripts per million) > 0.1. To download the list of genes and corresponding TCGA sample identifiers (across all 33 cancer types) with that gene expressed at TPM > 0.1, run the following: 

  ```bash
  EXPR_FILE="TCGA_GRCh38_expressed-genes_TPM.tsv.gz"
  wget http://compbio.cs.princeton.edu/pertinint/${EXPR_FILE}
  ```
  
  Of course, if you are looking at a different set of samples and do not have any expression information (or if you would prefer not to limit somatic mutations by expression), simply set the `--no_expression` tag when running PertInInt.

* **Known Driver Genes.** Once PertInInt is finished running, we include an option to annotate the resulting gene list with whether each gene is a known cancer driver gene. You can turn this option off with the `--no_driver_id` tag when running. To download the lists of known driver genes that we use for validation in our paper, run the following:

  ```bash
  DRIVER_FILE="GRCh38_known_drivers.tsv.gz"
  wget http://compbio.cs.princeton.edu/pertinint/${DRIVER_FILE}
  ```

  Note that you can customize this file to annotate your resulting gene list however you want. 

### 2: Run PertInInt

* PertInInt parses the input .maf file and stores mutation information in memory; you should allot enough RAM to the program to store this file. There are no further machine nor processor requirements. To run PertInInt: 

  ```bash
  python PertInInt.py --track_path track_weights/
                      --maf_file mafs/TCGA.Aggregate.muse.aggregated.somatic.maf
                      --expr_file TCGA_GRCh38_expressed-genes_TPM.tsv.gz
                      --out_file output/TCGA.Aggregate.muse.aggregated.somatic-PertInInt_output.tsv
  ```

### 3: Parsing PertInInt output

* By default, PertInInt returns a list of genes ordered in descending order by *Z*-score. Each gene is associated with a comma-delimited list of functional regions and their associated individual (positive) *Z*-scores. We provide a script for your convenience to parse this output and create a tab-delimited table highlighting which track types (e.g., specific interaction sites or domains) "light up" for particular genes to aid in downstream functional analyses. 

  ```bash
  python highlight_mechanism.py --pertinint_results output/TCGA.Aggregate.muse.aggregated.somatic-PertInInt_output.tsv
                                --out_file output/TCGA.Aggregate.muse.aggregated.somatic-PertInInt_mechanisms-output.tsv
  ```

### 4: Source code for customizing PertInInt tracks

Our analytical method can easily be extended to include other track types. The full source code as well as instructions for preprocessing tracks and precomputing the expectation, variance, and covariances needed to compute *Z*-scores can be found on the [PertInInt "internal" wiki](https://github.com/Singh-Lab/PertInInt-internal/wiki). 
