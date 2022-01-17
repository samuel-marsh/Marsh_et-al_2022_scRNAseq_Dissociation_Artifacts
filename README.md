# Marsh_et-al_2022_scRNAseq_Dissociation_Artifacts  
#### Code to reproduce analysis objects for the data contained in:  
[**Dissection of artifactual and confounding glial signatures by single cell sequencing of mouse and human brain**](LINK_TBD)  
Samuel E. Marsh<sup>1,\* </sup>, Alec J. Walker, Tushar Kamath<sup>1</sup>, Lasse Dissing-Olesen, Timothy R. Hammond<sup>2</sup>, T. Yvanka de Soysa, Adam M.H. Young, Sarah Murphy, Abdulraouf Abdulraouf, Naeem Nadaf, Connor Dufort, Alicia C. Walker, Liliana E. Lucca, Velina Kozareva<sup>2</sup>, Charles Vanderburg, Soyon Hong, Harry Bulstrode, Peter J. Hutchinson, Daniel J. Gaffney, David A. Hafler, Robin J.M. Franklin, Evan Z. Macosko, & Beth Stevens

<sup><sup>1</sup>Performed analysis</sup>   
<sup><sup>2</sup>Assisted analysis</sup>  
<sup><sup>\*</sup>Analysis lead (contact: samuel.marsh@childrens.harvard.edu)</sup>  

A copy of the code/repository which contained analyses from preprint can be downloaded [here](https://github.com/samuel-marsh/Marsh_et-al_2022_scRNAseq_Dissociation_Artifacts/tree/master/12_Code%20for%20Preprint%20Analyses).

## Code
Included is the code necessary to replicate the Seurat or LIGER (or both) objects used for analysis and plotting.
- Each R file specifies version of Seurat/LIGER used for analysis/object creation.
    - Some analyses were performed across multiple versions of Seurat (V2 > V3).  In this scenario objects were updated to V3 using `UpdateSeuratObject`
    - Scripts specify point of upgrade to V3 in regard to analysis or object modification.
    - Seurat V2.3.4 source package can be downloaded here from [CRAN Archive](https://cran.r-project.org/src/contrib/Archive/Seurat/) and installed from local source.
    - Seurat V3.2+ was released near the end of analysis.  To maintain consistency, Seurat V3.1.5 was downloaded from [CRAN Archive](https://cran.r-project.org/src/contrib/Archive/Seurat/) and installed from local source when switching between V2 and V3 was necessary.  

- Where possible date of analysis performed prior to is specified.  To replicate analyses performed on specific date the following actions are recommended or described in code:
  - Use of contained environment using [packrat](https://cran.r-project.org/web/packages/packrat/index.html) or [renv](https://cran.r-project.org/web/packages/renv/index.html) packages. Followed by date-specific version installation of CRAN packages using [versions](https://cran.r-project.org/web/packages/versions/index.html) package.
  - Archived source versions of specific packages may also be needed depending on version of R and can be downloaded from CRAN archives and installed from local source.

- LIGER analyses were performed using the in development ["online"](https://github.com/MacoskoLab/liger/tree/online) branch, updating throughout analysis to accommodate bug fixes.  
  - LIGER analyses also utilize multiple versions of Seurat as specified in code for some of the following situations:
    - Seurat V3 used used for data import, QC filtering (genes, UMIs, % mito), and majority of plotting.
    - Seurat V2 was used during LIGER analysis workflow to accommodate use of now deprecated [`clusterLouvainJaccard` function](https://github.com/samuel-marsh/Marsh_et-al_2020_scRNAseq_Dissociation_Artifacts/tree/master/11_Misc) which relied on Seurat V2 object structure.
    - Conversion between Seurat and LIGER objects was performed using built in LIGER functions `seuratToLiger` and `ligerToSeurat`.

## Data  
### Original Data
The data in this project can be broadly divided into 2 categories (7 sub-projects).  Please see [SI Table 1 & 2](ADD_NEW_LINK) (SI Table 1: Mouse Experiments 1-4) and (SI Table 2; Human Experiments 4-7 & Human Literature Reanalysis) for breakdown by sample and more information.

***A brief overview with links to the raw data (fastqs) and processed data (Cell Ranger `count` Gene Expression Matrices) see table below***
| Experiment | Species | Seq Used | Description | Raw/Count Data |
| :-----: | :-----: | :------: | :------------: | :---------: |
| Exp. 1 | Mouse | scRNA-seq (10X 3' V2) | scRNA-seq of microglia with 4 different dissociation protocols | [GSE152183](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152183) |
| Exp. 2 | Mouse | scRNA-seq (10X 3' V2) | scRNA-seq of all CNS cells with or without inhibitors | [GSE152182](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152182) |
| Exp. 3 | Mouse | scRNA-seq (10X 3' V2) | scRNA-seq of microglia (tail vein PBS injection) | [GSE152210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152210)|
| Exp. 4 | Mouse | scRNA-seq (10X 3' V3.0 & V3.1) | scRNA-seq of microglia w or w/o Inhibitors (10X Version Analysis) | [*GSE188441*](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188441)|
| |  |  | | |
| Exp. 5 | Human | snRNA-seq (10X 3' V3.0) | snRNA-seq of post-mortem brain tissue | [GSE157760](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157760)|
| Exp. 6 | Human | snRNA-seq (10X 3' V3.0) | snRNA-seq of surgically resected brain tissue with or without freezing time delay | [*in-progress*](EGAXXXXXXX)|
| Exp. 7 | Human | scRNA-seq (10X 5' V1) | scRNA-seq| [*phs002222.v2.p1*](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002222.v1.p1)|

### Processed Data
All proceesed data files represent the output from Cell Ranger `count`.  Files provided are the "*filtered_feature_bc_matrix*" (i.e. only containing the barcodes that Cell Ranger called as cells during preprocessing). Information on Cell Ranger version and Genome/Annotation for each experiment can be found in [SI Table 1 & 2](UPDATE_LINK) as well as individual repository meta data.  

**Experiments 1-4, 5 (NCBI GEO)**  
There are 3 processed data files per library:
  1. GSM\*\_*Sample-Name*_barcodes.tsv.gz: corresponds to the cell barcodes (i.e. column names).
  2. GSM\*\_*Sample-Name*_features.tsv.gz: corresponds to the gene identifiers (i.e. row names).
  3. GSM\*\_*Sample-Name*_matrix.mtx.gz: expression matrix in sparse format.

### Raw fastq Files
All raw data fastq/BAM files can be downloaded from SRA linked from NCBI GEO records, or from EGA/dbGaP records.

### Literature Reanalysis
Reanalyzed data from literature is summarized detailed in table below.
| Dataset | Species | Seq Used | Raw/Count Data | Publication |
| :-----: | :-----: | :------: | :------------: | :---------: |
| Mathys | Mouse | scRNAseq (Smart-seq2) | [GEO103334](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103334) <br> & Authors<sup>a</sup> | [Mathys et al., 2017 <br> (Cell Reports)](https://doi.org/10.1016/j.celrep.2017.09.039) |
| Plemel | Mouse | scRNAseq (10X 3' V2) | [GSE115803](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115803) | [Plemel et al., 2020 <br> (Science Advances)](https://doi.org/10.1126/sciadv.aay6324) |
| Zywitza | Mouse | scRNAseq (Drop-Seq) | [GSE111527](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111527) | [Zywitza et al., 2018 <br> (Cell Reports)](https://doi.org/10.1016/j.celrep.2018.11.003) |
| Mizrak | Mouse | scRNAseq (Microwell Seq) | [GSE109447](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109447) | [Mizrak et al., 2019 <br> (Cell Reports)](https://doi.org/10.1016/j.celrep.2018.12.044) |
| Zeisel | Mouse | scRNAseq (10X 3' V1 & V2) | [mousebrain.org](mousebrain.org) | [Zeisel et al., 2018 <br> (Cell)](https://doi.org/10.1016/j.cell.2018.06.021) |
| Hammond | Mouse | scRNAseq (10X 3' V1 & V2) | [GSE121654](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121654) | [Hammond et al., 2019 <br> (Immunity)](https://doi.org/10.1016/j.immuni.2018.11.004) |
| Keren-Shaul<sup>b</sup> | Mouse | MARS-Seq | [GSE98969](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98969) <br> | [Keren-Shaul et al., 2017 <br> (Cell)](https://doi.org/10.1016/j.cell.2017.05.018) |
| Pasciuto | Mouse | scRNAseq (10X 3' V2) | [GSE144038](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144038) <br> & [Mendeley Data](https://data.mendeley.com/datasets/hsmzw47kbg/3) | [Pasciuto et al., 2020 <br> (Cell)](https://doi.org/10.1016/j.cell.2020.06.026) |
| Crinier | Mouse | scRNAseq (10X 3' V2) | [GSE119562](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119562) | [Crinier et al., 2018 <br> (Cell)](https://doi.org/10.1016/j.immuni.2018.09.009) |
| Pasciuto | Human | scRNAseq (10X 3' V2) | [GSE146165](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146165) <br> & [Mendeley Data](https://data.mendeley.com/datasets/hsmzw47kbg/3) | [Pasciuto et al., 2020 <br> (Cell)](https://doi.org/10.1016/j.cell.2020.06.026) |
| Zhou | Human | snRNAseq (10X 5' V1) | [syn21670836](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage?Study=syn21670836) | [Zhou et al., 2020 <br> (Nature Medicine)](https://doi.org/10.1038/s41591-019-0695-9) |
| Morabito<sup>i</sup> | Human | snRNAseq (10X 3' V3.0) | [syn18915937](https://www.synapse.org/#!Synapse:syn18915937/wiki/592740) | [Morabito et al., 2020 <br> (Human Molecular Genetics)](https://doi.org/10.1093/hmg/ddaa182) |
| Leng & Li | Human | snRNAseq (10X 3' V2) | [syn21788402](https://www.synapse.org/#!Synapse:syn21788402/wiki/601825)<sup>c</sup> <br> & [GSE147528](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147528) | [Leng & Li et al., 2021 <br> (Nature Neuroscience)](https://doi.org/10.1038/s41593-020-00764-7) |

<sup><sup>a</sup>FPKM data and raw fastq files are available via GEO.  Raw count matrix was obtained via personal communication with authors.</sup>  
<sup><sup>b</sup>Only a specific subset of samples were used in reanalysis.  See reanalysis code for more information.</sup>    
<sup><sup>c</sup>Data on synapse are post-QC and were used for re-analysis.  GEO records contain the all barcodes (unfiltered) HDF5 cellranger output files and fastqs.</sup>  
<sup><sup>i</sup>Reanalysis of Morabito et al., was also used for calculation of cell type proportions in [Liddelow, Marsh, & Stevens et al., 2020 (Trends in Immunology)](https://doi.org/10.1016/j.it.2020.07.006)</sup>

### Human Data Reanalysis Meta Data
Meta data for human data was assembled from published SI Tables, public data on synapse, or restricted access data on synapse
  - Compiled publicly available meta data variables for each human dataset can be found in [SI Table 2](UPDATE LINK).
  - "DUC" in the table indicates data available from synapse following submission and approval of Data Use Certificate.

## Acknowledgements:
This study was supported by funding from Cure Alzheimer's Fund (B.S.).  Special thanks to authors Tushar Kamath, Tim Hammond, Alec Walker, Lasse-Dissing-Olesen, Velina Kozareva, Evan Macosko, as well other members of Stevens and Macosko labs for helpful discussions and assistance during the analysis phase of this project.  

**Data Acknowledgements:**  
The analysis and results published here from Zhou et al., 2020 in whole or in part are based on data obtained from the [AMP-AD Knowledge Portal](https://adknowledgeportal.synapse.org/). Samples for this study were provided by the Rush Alzheimer’s Disease Center, Rush University Medical Center, Chicago. Data collection was supported through funding by NIA grants P30AG10161, R01AG15819, R01AG17917, R01AG30146, R01AG36836, U01AG32984, U01AG46152, the Illinois Department of Public Health, and the Translational Genomics Research Institute. Raw data used in analysis here are available from AMP-AD/Synapse database through links provided in table above.  Additional ROSMAP data can be requested at [https://www.radc.rush.edu](https://www.radc.rush.edu).

The analysis and results published here for Morabito et al., 2020 are based on reanalysis of study data downloaded from Synapse as provided by Dr. Vivek Swarup, Institute for Memory Impairments and Neurological Disorders, University of California, Irvine.  Data collection was supported through funding UCI Startup funds and American Federation of Aging Research.  Raw data used in analysis here are available from the Synapse database through link provided in table above.
