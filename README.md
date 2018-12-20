# PertInInt

<div align="center"><img src="http://compbio.cs.princeton.edu/pertinint/PertInInt_workflow.png" alt="pipeline figure" title="PertInInt workflow" width="80%" /></div><br />

The goal of our analytical PertInInt method is to rapidly uncover proteins with significant enrichments of somatic **Pert**urbations **In** **Int**eraction and other functional sites. If you use data or scripts from this repository, please cite:

> Kobren, S.N., Chazelle, B. and Singh, M. (2019). "An integrative approach to identify preferentially altered interactions in human cancers." *Manuscript in preparation.*

### 1: Downloading preliminary data

* PertInInt uses *tracks* to model different functional regions of a protein. These tracks measure interaction interfaces, functional protein domains and conserved protein positions. To download the set of precomputed tracks used by PertInt (which will be unzipped into a folder called `track_weights`), run the following: 

  ```bash
  PERTININT_TRACKS="PertInInt-tracks_v0.tar.gz"
  wget http://compbio.cs.princeton.edu/pertinint/${PERTININT_TRACKS}
  tar -xvzf ${PERTININT_TRACKS}
  ```

* PertInInt runs on any input .maf file. We tested PertInInt using a .maf file containing somatic mutations from all 33 TCGA cancer types, obtained from [NCI's Genomic Data Commons](https://gdc.cancer.gov) on December 6, 2018. To download this aggregated .maf file to test on, run the following: 

  ```bash
  if [ ! -d mafs ]; then mkdir mafs; fi
  AGGREGATE_CANCER="TCGA.Aggregate.muse.aggregated.somatic.maf.gz"
  wget http://compbio.cs.princeton.edu/pertinint/${AGGREGATE_CANCER} -O mafs/${AGGREGATE_CANCER}
  gzip -d mafs/${AGGREGATE_CANCER}
  ```

### 2: Run PertInInt

* PertInInt parses the input .maf file and stores it in memory; you should allot as much memory as is required here to run PertInInt. This program can be run on any machine and does not have specific processor requirements. To run PertInInt: 

  ```bash
  python PertInInt.py --track_path track_weights/
                      --input_maf mafs/TCGA.Aggregate.muse.aggregated.somatic.maf
                      --output output/TCGA.Aggregate.muse.aggregated.somatic-PertInInt_output.tsv
  ```

### 3: Parsing PertInInt output

* By default, PertInInt returns a list of genes ordered in descending order by *Z*-score. Each gene is associated with a comma-delimited list of functional regions and their associated individual (positive) *Z*-scores. We provide a script for your convenience to parse this output and create a mechanism plot highlighting which track types "light up" for particular genes to aid in functional analyses. 

