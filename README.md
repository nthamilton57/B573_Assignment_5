# B573 Assignment 5: Advanced R

## Programming Language: R

## Date: November 30, 2023

## Description:

This repository is designed to servee as an extension to learning R in the realm of bioinformatic analyses. Through the assignment, package creation, the package creation workflow, and vignette creation were used. Additionally, bioinformatic analyses like differential gene expression were explored to see how real world datasets can be analyzed. The necessary prerequisite steps were explored as well, such as dimensionality reduction. Finally, extensions of the expression analysis were explored.


### Required Files


<code>Assignment_5.R</code>  retrieves gene data from GEO, prepares it, and does differential gene expression analysis

<code>GeoDE Vignette.rmd</code>  creates a vignette displaying the functionality of the GeoDE package


### Required Packages


<code>Assignment_5.R</code>  forcats, stringr, ggplot2, ggrepel, readr, tidyr, survminer, dplyr, BiocManager: GEOquery, limma, pheatmap, org.Hs.eg.db, ReactomePA, DOSE

<code>GeoDE Vignette.rmd</code>  GeoDE (downloaded from this GitHub)


### Execution

<ins>Steps to Run the R File</ins>

1. Clone repository

2. Install required dependencies

3. Change working directory to location of files

    ```
    cd path/to/files
    ```

4. Run the Assignment_5 script using the following commands: 

    ```
    Rscript assignment_5.R
    ```

5. The script will output results in the command line.

6. Using each of the functions in the <code> GeoDE Package</code>  file, run the <code>GeoDE Vignette</code> R markdown file.


### Output Files

No Output Files
