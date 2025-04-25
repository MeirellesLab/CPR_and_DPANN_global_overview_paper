# CPR and DPANN Global overview paper

This repository stores the data and scripts for all steps, analyses and figures of the paper by Ferreira et al. (2025) "Low-Abundance Microorganisms and Their Role in Driving Microbial Community Structure Dissimilarities Across Ecosystems".

Below is a breakdown of each step taken to obtain the results of the article:

### 1. Samples download

### 2. Kraken database assembly

### 3. Metagenome annotation

### 4. Analysis and figures

#### 4.1. Run summarizations and analysis

 - The script "general_prospect.R" contains all the code that generates the analyses and summarization for figures 1, 2 and 3, as well as for some supplementary tables and also for SIMPER, which is used in figure 4.
 - In order, he does the following:
    - Sample counting for Table S6;
    - General NMDS analysis and General PERMANOVA analysis (lifestyle and ecosystem)
        - For figure 1
    - SIMPER analysis (ecosystem and habitat)
        - For figure 4c
    - Assignment of microbial groups (Bonafide, CPR, DPANN), calculation of site-level abundance and richness, summarization of microgroup prevalence, prevalence summaries by lifestyle and ecosystem.
        - For figure 3
    - NMDS analysis by microbial group (Bonafide, CPR, DPANN) and PERMANOVA analysis by microbial group (Bonafide, CPR, DPANN)
        - For figure 1
    - GAM models for abundance and richness vs. latitude
        - For figure 2
    - PerANOVA analysis for abundance and richness by group
        - For figure 3

#### 4.2 Run scripts to generate figures
 - The script "figures_script.R" contains all the code to generate figure 1, 2 and 3 for the paper. It utilizes the RData generated for the analysis in the previous step and generate the figures.

 - The script "heatmap_jitter.R" contains all the code to generate figure 4 and figure S5. 

 _Observation to following work: I need to find the scripts that generate the other supplementary figures to describe here_

 
