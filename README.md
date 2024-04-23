# fMRI-data-analysis
From preprocessing to RSA
# MRI Data Analysis Scripts

This package contains several MATLAB scripts used for the processing and analysis of MRI data, specifically designed for functional magnetic resonance imaging (fMRI). Each script is tailored for specific steps in the MRI data processing pipeline:

## Scripts Overview

### `SRC_Preprocessed.m`
- **Purpose:** Preprocesses MRI data for further analysis.
- **Key Features:** Includes steps such as reading NIfTI files, performing slice timing correction, motion correction, and normalization of images to a standard brain template using SPM (Statistical Parametric Mapping).

### `SRC_FstLeve.m`
- **Purpose:** First-level statistical analysis in the context of fMRI.
- **Key Features:** Loads and processes condition and fixation data, handles the setup of fMRI analysis parameters in SPM, and manages file paths for input/output data.

### `SRC_RSA_ROIbased.m`
- **Purpose:** Executes ROI-based representational similarity analysis.
- **Key Features:** Uses predefined regions of interest to extract activity patterns and compute similarity matrices.

### `SRC_RSApermutation1st.m`
- **Purpose:** Manages the permutation tests for RSA based on ROIs, typically used to assess statistical significance of results.
- **Key Features:** Implements multiple rounds of permutation to create a distribution of chance-level results for comparison with actual data.

### `SRC_Searchlight_Fusiform_RSA.m`
- **Purpose:** Conducts searchlight analysis across the fusiform gyrus for detailed RSA.
- **Key Features:** Standard RSA without permutation tests.

### `SRC_Fusiform_SearchlightRSAperm.m`
- **Purpose:** Performs a searchlight RSA in the fusiform gyrus with permutation testing.
- **Key Features:** Executes fine-grained analysis by moving a searchlight across the brain to identify areas where activity patterns correlate with experimental conditions.

  ### `SRC_WhoBraRSA.m`
- **Purpose:** This script handles the whole-brain RSA (representational similarity analysis), comparing neural patterns across different conditions within subjects.
- **Key Features:** Computes representational dissimilarity matrices.

## General Information
- **Dependencies:** All scripts require MATLAB and the SPM12 toolbox.
- **Installation:** Ensure that all file paths are correctly set in each script to match the local setup.
- **Usage:** Each script is executed individually, depending on the specific stage of analysis required.

Feel free to modify and adapt these scripts as needed for your specific research needs. For any issues or further information, contact the original author of the scripts.

----
reference
Hebart, M. N., Görgen, K., & Haynes, J.-D. (2015). The Decoding Toolbox (TDT): a versatile software package for multivariate analyses of functional imaging data [Methods]. Frontiers in Neuroinformatics, 8. https://doi.org/10.3389/fninf.2014.00088 
Liu, X., Wisniewski, D., Vermeylen, L., Palenciano, A., Liu, W., & Brysbaert, M. (2021). The Representations of Chinese Characters: Evidence from Sublexical Components. The Journal of Neuroscience, 42, JN-RM. https://doi.org/10.1523/JNEUROSCI.1057-21.2021 
Oganian, Y., Conrad, M., Aryani, A., Spalek, K., & Heekeren, H. R. (2015). Activation Patterns throughout the Word Processing Network of L1-dominant Bilinguals Reflect Language Similarity and Language Decisions. J Cogn Neurosci, 27(11), 2197-2214. https://doi.org/10.1162/jocn_a_00853 
Patel, T., Morales, M., Pickering, M. J., & Hoffman, P. (2023). A common neural code for meaning in discourse production and comprehension. NeuroImage, 279, 120295. https://doi.org/https://doi.org/10.1016/j.neuroimage.2023.120295 


