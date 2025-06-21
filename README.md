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
   3. 2a_pc_detection -- qualitatively determine the number of PCs to use for downstream analysis (at elbow)
   4. 2b_global_search -- (search only) quantitatively determine the number of features and clustering resolution for downstream analysis
   5. 2c_global_search_postanalysis -- determines the parameters for 2b
   6. 3a_clustering -- use the parameters above for clustering
   7. 3b_clustering_sensitivity_test -- test the stability and sensitivity of the clustering
   8. 3c_clustering_pvalue -- see if clustering if different from noise
   9. 4_figures_main -- figure plotting for the paper
   10. 5_figures_supp_donor -- supplementary figures comparing animal and species affects
3. Figures will be saved into /analysis_output folder
4. Optional scripts:
   - sample_umaps -- quickly test your clustering and umaps with this

## Contact
Scott Sun - shsun@mgh.harvard.edu
John Pezaris - pezaris.john@mgh.harvard.edu
