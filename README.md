# Bolaji_Luthi_RNAseq

Used for basic alignment and exploratory analysis of RNAseq data for:

**Cohesin forms fountains at active enhancers in C. elegans** *Bolaji N. Lüthi, Jennifer I. Semple, Anja Haemmerli, Saurabh Thapliyal, Kalyan Ghadage, Klement Stojanovski, Dario D’Asaro, Moushumi Das, Nick Gilbert, Dominique A. Glauser, Benjamin Towbin, Daniel Jost, Peter Meister bioRxiv 2023.07.14.549011; doi: https://doi.org/10.1101/2023.07.14.549011*

[now accepted in Nature Communications]

Based on pipeline created by Jenny Semple (SMC_RNAseq) and Todor Gitchev(CeFTALL).

Final figures used in paper are in https://github.com/CellFateNucOrg/Luthi_etal

## Installation and preparation

This only needs to be done once

### 0.1 Clone project:

From command line:

    git clone git@github.com:CellFateNucOrg/Bolaji_RNAseq.git

    # how to update the remote url when iusing access token:
    git remote set-url origin git@github.com:CellFateNucOrg/Bolaji_RNAseq.git
    # how to see what is current set:
    git config --get remote.origin.url
    # after propperly set accessible remote url, you should be able to do: git pull 

### 0.2 Python conda:

Python Miniconda installation if needed:  
```
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh
    ./Miniconda3-latest-Linux-x86_64.sh -b
    
    # to use conda with bash shell 
    conda init
```
Add a CONDA_ACTIVATE variable pointing to conda activate binary in your ~/.bashrc script, as it doesn't always work in SLURM environment:

```
echo "export CONDA_ACTIVATE=${HOME}/miniconda3/bin/activate" > ~/.bashrc
```

### 0.3. Install environment:

Installation will be done primarily with conda you need to make sure the base environment is up to date. 

run script: 
```
    ./install_env.sh
```

For R installation there is some variation between bioinformatics and ubelix clusters where you can use module load and Pertz cluster where you use a singularity image of R.
  

### 0.4. Downloads reference data and index genome/transcripts

Change location of where to put reference genome in _00_downloadAndIndexGenome.sh
run script:

    ./_00_downloadAndIndexGenome.sh

   
## Run RNA-seq analysis

### 1.1 Prepare your data

Place all desired .fastq.gz files in a folder and create a <fastqList>.csv with the column structure:

    #fastqFile,sampleID,strain
    
If you have PE sequencing, use the following headers:
    
    #fastqFile1,fastqFile2,sampleName,strain
    
Names of the fastqFiles (with full path) and the sampleName (sampleName must be unique to each row) are essential for mapping. The other fields are used by DESeq2 to create appropriate comparison groupings and can be changed according to your data.

### 1.2 Map RNA-seq reads

Open the _01_mapRNA.sh file and change the #SBATCH --array=1-24%5 line to reflect the number of samples in your fastqList.csv file (here 24 files which will be processed in batches of 5 so as not to overload the server).

Make sure fastqFileList variable is set to your fastqList.csv file name
Make sure the genomeVer and GENOME_DIR variables are correctly set (same as in the indexing script)

In SLURM environment use the following command:

```
sbatch _01_mapRNAreads.sh 
```    

### 1.3 Differential Expression (DESeq2) analysis

Differential expression and some QC plots are performed in R running the following script:

```
DESeq2_analysis_2_coh1.R
```

As QC, compared results with DESeq2 results from Das et al. where same RNAseq was processed together with cleavage of other SMC complexes. Data slightly different, but largely similar (correlation 0.97-99). See *./compareDatasets.R* script.

### 1.4 Further analysis

NOTE: Many libraries used in these scripts require installation in R as they were not included in the setup files.

**Enhancers:**

*./compareGenesWithEnhancers.R*


**Fountain features:**

*./compareAtFountains.R*


**GO and tissue enrichement:**

*./compareTissueAndGO.R*

*./allPlots_TEA.Rmd*

*./allPlots_WORMCAT.Rmd*


**Compare to single cell RNAseq:**

*./compareToScRNAseq.R*


**Gene set enrichment and orthologs in CdLS:**

*./compareGeneLists_diopt_WormToHuman.R*

*./compareGeneLists_ortholist_WormToHuman.R*

*./compareGeneLists_ortholist_HumanToWorm.R*

### 2. Differentail transcript expression (Sleuth)

*./Sleuth_analysis.R*

