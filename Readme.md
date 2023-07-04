A 4D Flow Processing Module Based on [Quantitative Velocity Tool (QVT)](https://github.com/uwmri/QVT)
=========
## Introduction ##
Built on QVT, some extra functionality has been added to extract time continuous boundary conditions which can be used to drive flow in cerebral hemodynamic models. 

Major changes from QVT include the development of an improved flow consistency metric which allows the automatic extraction of the ideal flow from each cerebral vessel. This was necessary as there are several "good" flow positions which increase the interobserver error when using this flow profile as input into patient specific modelling. 

We also include the ability to import DICOM and NIFTY files instead of the .dat required in original QVT.

Other functions convert the flow into periodic Fourier coefficients 

UPDATE more as this is developped.


## Installation ##
Requires MATLAB version > 2018
### Dependencies ###
**Required Matlab Add-Ons** \
Image Processing Toolbox (for medfilt3) \
Curve Fitting Toolbox (for csaps) \
Statistics and Machine Learning Toolbox (for kmeans)
## Usage
Run paramMap as per QVT, once the branch has loaded, click any point for the branches of interest and save to the generated csv. 

Run the separate function "flow_pipeline" which will let you point to where the subject data is saved, it will take the branch ID's, and then extract the quality function over the branch and output the highest quality flows automatically. The average flow (based on quality cutoff) is used for boundary conditions. Plots of flow and error are also included in a separate csv. 


### Citations ### 
If you are using this module, be sure to cite the original creators and publications:

- [Roberts GS, Hoffman CA, Rivera-Rivera LA, Berman SE, Eisenmenger LB, Wieben O. Automated hemodynamic assessment for cranial 4D flow MRI. Magn Reson Imaging. 2022 Dec 26:S0730-725X(22)00231-4. doi: 10.1016/j.mri.2022.12.016. Epub ahead of print. PMID: 36581214](https://pubmed.ncbi.nlm.nih.gov/36581214/)

- [Schrauben E, Wahlin A, Ambarki K, Spaak E, Malm J, Wieben O, Eklund A. (2015). Fast 4D flow MRI intracranial segmentation and quantification in tortuous arteries. J Magn Reson Imaging, 42(5), 1458-1464. doi:10.1002/jmri.24900](https://pubmed.ncbi.nlm.nih.gov/25847621/)
