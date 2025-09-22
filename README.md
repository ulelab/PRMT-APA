
**PRMT-APA README**

This repository contains the scripts and data files necessary to generate the figures from "PRMT activity promotes global 3' UTR shortening in proliferating cells" (Griffith et al., 2025, doi: https://doi.org/10.1101/2025.03.06.641848)



**RAW DATA**

The data generated in this study has been deposited at EBI ArrayExpress under the following accession numbers: 
- “E-MTAB-14928” (96hrs ADMAi vs SDMAi vs SDMAi+ADMAi)
- “E-MTAB-14927” (DMAi timecourse)
- “E-MTAB-14932” (DMAi across panel of 10 cancer lines)
- “E-MTAB-14990” (murine T cells)
- “E-MTAB-14992” (patient-derived organoids)
- “E-MTAB-14930” (siCFIm25 +/- DMAi)

The fastq files present in each of these datasets were processed using the nf-quantseq pipeline (v.0.1.0; DOI: 10.5281/zenodo.14417140), which generates Count bed files and a poly(A) site atlas - all of which are located in the Data/ directories in this repository, organised by Figure. These files were used as input for the DRIMseq.R scripts which were the starting point for all APA analyses performed in the manuscript.



**REPRODUCING ANALYSES FROM MANUSCRIPT**

For the more complex scripts located in this repository, script-specific READMEs can be found in the same directories as those scripts. Files required for each script are located in the Data/ directory also found in this repository, organised by figure. Following the initial DRIMSeq analysis, APA events were next classified using the classification_of_APA_events.R script, before all other dataset-specific secondary analyses were peformed.



**CONTACT**

Contact llywelyn.griffith@kcl.ac.uk if you have any questions.
