# PATTY: a computational method for correcting open chromatin bias in bulk and single-cell CUT&Tag data 

Precise profiling of epigenomes is essential for better understanding chromatin biology and gene regulation. Cleavage Under Targets & Tagmentation (CUT&Tag) is an easy and low-cost epigenomic profiling technique that can be performed on a low number of cells and at the single-cell level. With its growing adoption, CUT&Tag datasets spanning diverse biological systems are rapidly accumulating in the field. CUT&Tag assays use the hyperactive transposase Tn5 for DNA tagmentation. Tn5’s preference toward accessible chromatin alters CUT&Tag sequence read distributions in the genome and introduces open chromatin bias that can confound downstream analysis, an issue more substantial in sparse single-cell data. We show that open chromatin bias extensively exists in published CUT&Tag datasets, including those generated with recently optimized high-salt protocols. To address this challange, we present PATTY (Propensity Analyzer for Tn5 Transposase Yielded bias), a comprehensive computational method that corrects open chromatin bias in CUT&Tag data by leveraging accompanying ATAC-seq. By integrating transcriptomic and epigenomic data using machine learning and integrative modeling, we demonstrate that PATTY enables accurate and robust detection of occupancy sites for both active and repressive histone modifications, including H3K27ac, H3K27me3, and H3K9me3, with experimental validation. We further develop a single-cell CUT&Tag analysis framework built on PATTY and show improved cell clustering when using bias-corrected single-cell CUT&Tag data compared to using uncorrected data. Beyond CUT&Tag, PATTY sets a foundation for further development of bias correction methods for improving data analysis for all Tn5-based high-throughput assays.

Our manuscript is now published on [Nature Communications](https://www.nature.com/articles/s41467-026-73599-8)<br>

## 0. Introduction to the PATTY Package

**PATTY** is a computational tool designed to correct open chromatin bias in CUT&Tag data at both **bulk** and **single-cell** levels. Current version of PATTY support open chromatin bias correction for H3K27me3, H3K27ac, and H3K9me3. It leverages a pre-trained logistic regression model, built using CUT&Tag data in the K562 cell line, to correct bias for specific histone modifications.

- **Bulk mode:** PATTY applies the correction model to genome-wide 200bp tiling bins, and generates a bias-corrected score on each candidate bin. 

- **Single-cell mode:** PATTY performs bias correction at the individual cell level, producing a 200bp-bin by cell matrix of corrected signals. It then supports downstream cell clustering analysis using the bias-corrected data to improve biological interpretability and resolution.


- Changelog<br>
v1.0 PATTY for biorxiv manuscript and initial submission<br> 
v1.1 Improve the installation steps. Designed for paper revision.

## 1. Installation
- Package requirements<br>
PATTY requires Linux or MacOS as OS.<br>
PATTY requires [Python](https://www.python.org) 3.6+ and [Rscript](https://www.r-project.org) v3+ to run.<br>
PATTY requires Python packages [scipy](https://scipy.org), [numpy](https://numpy.org), [scikit-learn](https://scikit-learn.org/stable/), and [joblib](https://joblib.readthedocs.io/) pre-installed. <br>
PATTY requires Bedtools (Quinlan et al., Bioinformatics, 2010) and UCSC tools (Kuhn et al., Brief Bioinform. 2013) pre-installed. <br>
PATTY sc mode requires Rpackage [ArchR](https://www.archrproject.com/), [pandas](https://pandas.pydata.org/), [tabix](http://www.htslib.org/doc/tabix.html), and [bgzip](http://www.htslib.org/doc/bgzip.html) pre-installed.

- Genome-wide mappable region annotation<br>
The genome-wide annotation file for [hg38](https://www.dropbox.com/scl/fi/k6iwvr0rh7ozycxv8u8fh/hg38_mappableBin.bed.gz?rlkey=5ofxbeankzvjlb608p2o7y3bn&st=x6y57wcx&dl=0) and [mm10](https://www.dropbox.com/scl/fi/1fa3ji3pmurreb05lopos/mm10_mappableBin.bed.gz?rlkey=250svcejszkf80zmjel8yaz5m&st=jvzu39ld&dl=0) genome can be downloaded here and input when running PATTY.

- Download PATTY <br>
\# Download PATTY from this github repo
```sh
$ git clone https://github.com/Tarela/PATTY.git
```

- Install PATTY <br>
# for root user
```sh
$ cd PATTY
$ pip install .  
```
\# if you are not the root user, you can install PATTY at a specific location where you have write permission
```sh
$ pip install . --prefix /home/PATTY  # Here you can replace “/home/PATTY” with any location 
$ export PATH=/home/PATTY/bin:$PATH    # setup PATH for the software
$ export PYTHONPATH=/home/PATTY/lib/python3.6/site-packages:$PYTHONPATH    # setup PYTHONPATH for module import
```
\# To check the PATTY package, just type:
```sh
$ PATTY --help  # If you see the help manual, you have successfully installed PATTY
```

\# NOTE: 
- To install PATTY on MacOS, the users need to download and install Command Line Tools beforehand

## 2. Run PATTY (usage)
#### Essential parameters
To run PATTY with the default parameters, you can set the following parameters:
-   -m MODE, --mode=MODE
Mode of PATTY, choose from bulk or sc(single-cell)
-   -c CUTTAG, --cuttag=CUTTAG
CUTTAG Input fragments file in (paired/single end) bed format for CUT&Tag data, with .bed extension (or .bed.gz for comparessed file). For sc mode, the 4th(name) column of the bed file represents the name/barcode of the corresponding individual cell
-   -a ATAC, --atac=ATAC
ATAC Input fragments file in bed format for ATAC-seq data, with .bed extension(or .bed.gz for comparessed file). The ATAC-seq fragments were used as bulk data for both sc and bulk modes (only chrm,start,end 3columns are required)
-   -f FACTOR, --factor=FACTOR
FACTOR Factor type of the CUT&Tag data. Currently PATTY support H3K27me3 (default), H3K27ac, and H3K9me3
-   -g GENOME, --genome=GENOME
genome version of the input data, choose from hg38 (default) and mm10
-   -o OUTNAME, --outname=OUTNAME
Name of output results

Example of running PATTY with default parameters (test data downloadable in :

\# bulk mode (mappable bin file can be downloaded in section1)
```sh
$ PATTY -m bulk -c ${path}/testbulk_H3K27me3.bed.gz -a ${path}/testbulk_ATAC.bed.gz -f H3K27me3 -o testbulk -g hg38 -b hg38_mappableBin.bed.gz 
```

\# sc mode (note that sc mode require bgzip compressed bed reads file for scCUT&Tag reads) 
```sh
$ PATTY -m sc -c ${path}/testsc_H3K27me3.bed.gz -a ${path}/testsc_ATAC.bed.gz -f H3K27me3 -o testsc -g hg38 -b hg38_mappableBin.bed.gz 
```


## 4. Pre-processing Steps for Generating the Input Fragments File

PATTY takes aligned fragment files in **BED format** as input(or .bed.gz for gzip comparessed file). Users may apply any preferred pre-processing pipeline to generate these files. We recommend retaining only **high-quality reads** with **MAPQ > 30** to ensure accurate bias correction. Note that PATTY takes original fragments bed files as input (e.g., transformed directly from aligned BAM files, or 10x cell ranger outputed fragments.tsv file for sc data). Please don't do any customized extension or shifting. 

### Default Input Format

The expected BED format varies depending on data type:

#### • Bulk CUT&Tag (Single-End)
```
chr1    10500   10646   .   .   +
chr2    20840   20986   .   .   -
```
> The 4-5th column represents an optional placeholder.

#### • Bulk CUT&Tag (Paired-End)
```
chr1    10500   10646
chr2    20840   21000
```

#### • Single-Cell CUT&Tag
```
chr1    10500   10646   CellA
chr2    20840   21000   CellB
```
> The 4th column must contain the cell barcode or cell name (like AATAACTACGCC-1).

## 5. Output Files

### Bulk Mode Outputs

1. `NAME_PATTYscore.bw`  
   A 200bp-resolution genome-wide track in **bigWig** format containing the PATTY scores for each candidate bin.  
   - Scores range from 0 to 1. Higher scores indicate higher confidence of true histone mark occupancy, while lower scores reflect likely false-positive or background signals due to open chromatin bias.

### Single-Cell Mode Outputs

1. `NAME_binXcell.txt.gz`  
   A **bin-by-cell PATTY score matrix** generated from single-cell CUT&Tag analysis.  
   - Rows: 200bp bins  
   - Columns: individual cells
   - Values: Similar PATTY score like in Bulk mode but for each individual cell  


## 6. Testing data and example of output files
We provided the test data for users to test PATTY. The sc/bulk output can also be generated with the command lines in Section 2 using the testing data as input. Click the file names to download. 
- testing data for **bulk** mode:
   - H3K27me3 [`Dropbox`](https://www.dropbox.com/scl/fi/7g46rtmeg9eihj0jbil08/testbulk_H3K27me3.bed.gz?rlkey=vgedovvqti7blf5d8bpriav7v&st=ozlmcxjm&dl=0)
   - ATAC [`Dropbox`](https://www.dropbox.com/scl/fi/0ls8fb7qw5v91vmjeff9h/testbulk_ATAC.bed.gz?rlkey=tkw4m8xzpe58m99foersrcpxq&st=dvsycypq&dl=0)
- testing data for **sc** mode:
   - H3K27me3 [`Dropbox`](https://www.dropbox.com/scl/fi/2pid9491b99jgzp784mzo/testsc_H3K27me3.bed.gz?rlkey=06orl5ot84q6r0p5zfxzuyzwm&st=8xfvowsj&dl=0)
   - ATAC [`Dropbox`](https://www.dropbox.com/scl/fi/osv616pafwoycjk5cfkx8/testsc_ATAC.bed.gz?rlkey=xqnyzdy65v9asnlrozh0vjwy6&st=rko7bzte&dl=0)


## 7. Other parameters in the PATTY pipeline
You can also set the following parameters for more accurate bias estimation and correction:
- -\-binMinReads=BINMINREADS  
[optional] Bins with < 5(default) reads covered will be discarded in the analysis. For sc mode, bins with a total of < 5 (default) reads across all high-quality cells will be discarded. set 0 to turn off this parameter. 
- -\-cellnames=CELLNAMES  
[optional] Single column plain text file for name list of used individual cells, each line contain the name of the individual cell. This parameter is only used for sc mode. 
- -\-readCutoff=READCUTOFF  
[sc optional] Reads number cutoff for high-quality cells. Cells with < 10000(default) reads will be discarded in the analysis. Users can change this parameter for samples with low sequencing depth to include more cells in the analysis. Setting a lower number for this parameter may decrease the accuracy of clustering results due to the low-quality cells. 
- -\-clusterMethod=CLUSTERMETHOD  
[sc optional] Method used for single-cell clustering analysis. The default is K-means (PCA dim reduction + K-means clustering). Optional choices (Seurat and scran) require related packages installed (described in section x)
- -\-clusterNum=CLUSTERNUM  
[sc optional] Number of clusters specified for K-means clustering and only used for the PCAkm (setting by --clusterMethod) method. The default is 7. 
- -\-UMAP  
[sc optional] Turn on this parameter to generate a UMAP plot for the clustering results.
- -\-overwrite  
[optional] Force overwrite; setting this parameter will remove the existing result! PATTY will terminate if there is a folder with the same name as -o in the working directory. Set this parameter to force PATTY to run. 
- -\-keeptmp  
[optional] Whether or not to keep the intermediate results (tmpResults/)

## 8. Reproduce figure results from the PATTY manuscript
Users can reproduce the bias correction results from the manuscript (Figure 4A, G, H, H3K27me3 CUT&Tag rep1) by running PATTY with the following command line:
```sh
$ PATTY -m bulk -c ${path}/H3K27me3_CUTTag_rep1.bed.gz -a ${path}/ATAC.bed.gz -f H3K27me3 -o testbulk -g hg38 -b hg38_mappableBin.bed.gz
```
Download input [CUT&Tag](https://www.dropbox.com/scl/fi/nyonvtpe8pdhm607vd151/H3K27me3_CUTTag_rep1.bed.gz?rlkey=liligwf0mt3mnyj678zcaloub&st=szd8yt24&dl=0) and [ATAC](https://www.dropbox.com/scl/fi/gxqxljiv6mr68l92e4hsm/ATAC.bed.gz?rlkey=n05eoajmo54a9lovz0bjba3rj&st=8nq2wcoa&dl=0) data, and example [output](https://www.dropbox.com/scl/fi/cavku8dzgydwj1t26wt9w/testbulk_correctSig.bw?rlkey=buvvkarm1ongmiwzny9z1d0ht&st=zsixclei&dl=0) here.

### PATTY Manuscript Source Data and Figure Generation Code
Users can reproduce other results and figures with the scripts in this [Reproduction Instructions](https://github.com/zang-lab/PATTY_figure_scripts).
