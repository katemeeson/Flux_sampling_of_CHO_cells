# **Flux sampling of CHO cells**
Flux Sampling Suggests Metabolic Signatures of High Antibody-Producing CHO Cells
## Link to publication in Biotechnology and Bioengineering: https://doi.org/10.1002/bit.28982
![image](https://github.com/user-attachments/assets/e7501eb9-ebda-43a7-8748-84dc5cfa900d)

### **Corresponding Author:** Kate Meeson: kate.meeson@manchester.ac.uk
## **Table of contents** ##
1. src
- Contains code to integrate transcriptomics data with GEM and perform flux sampling.
2. figures
- Contains code to generate figures from publication.
3. README.md
4. requirements.txt
5. RStudio_supplementary.txt
## **Repository layout** ##
```
├── src/                                                 # Folder containing source code for flux sampling and CBM construction
|  ├── 2_runFluxSampling_CSF.py                          # Integrate transcriptomics with GEM and run flux sampling
|  ├── function_mapping_transcriptome_data_JW.py         # Function to integrate transcriptomics with GEM
|  ├── iCHOv1_221-107_producing.ipynb                    # Modify iCHO2441 for FUJIFILM-specific protein product and validate model
├── figures/                                             # Folder with code to generate figures from publication (https://doi.org/10.1002/bit.28982)
│   ├── FluxSampling_analysis_and_visualisation.ipynb    # Study flux sampling solutions and identify high IgG-producing solutions. Generate Fig 4b and 4c
│   ├── FluxSampling_highIgG_reactions_subsystemsVis.R   # Subsystem enrichment analysis and generate Fig 5b
│   ├── DeSeq2_analysis.R                                # Generate R version of Fig 2a, 2b and 2c
│   ├── transcriptomics_figs.ipynb                       # Generate Python versions of Fig 2a, 2b and 2c
├── README.md
├── requirements.txt
├── RStudio_supplementary.txt                            # R-equivalent of requirements.txt file
```
## **Acknowledgements** ##
The authors acknowledge financial support from a Prosperity Partnership grant (EP/V038095/1) funded by EPSRC, BBSRC, and FUJIFILM Diosynth Biotechnologies. We also acknowledge everyone who has helped interpret experimental data for input datasets; those groups that developed the GEMs that are used here and the more general systems biology field for their support and feedback. We thank the Universities of Manchester, York and Edinburgh.
