[![DOI](https://img.shields.io/badge/DOI-10.1101/2024.11.14.623611-blue)](https://doi.org/10.1101/2024.11.14.623611)

# Description
R scripts to replicate the macaque LGN transcriptomics analysis results presented in "Single-cell transcriptomic analysis of macaque LGN neurons reveals novel subpopulations", [Sun et al., 2024](https://doi.org/10.1101/2024.11.14.623611).

## Updates
- 2025/6
   - Added global search for parameter optimization
   - Added analysis testing with ground truth data
   - Added sensitivity and stability testing

## Steps
1. Download macaque LGN data from the [Allen Brain Map](https://portal.brain-map.org/atlases-and-data/rnaseq/comparative-lgn) and place into /data folder
2. Run the R scripts in the following order
   1. 0_install_packages
   2. 1_data_collection
   3. 2_clustering
   4. 3_figures
3. Figures will be saved into /analysis_output folder

## Contact
Scott Sun - shsun@mgh.harvard.edu
John Pezaris - pezaris.john@mgh.harvard.edu
