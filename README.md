# FDCFIT: Fitting of flow duration curves using parametric formulations

## Description

The flow duration curve (FDC) is a signature catchment characteristic that depicts graphically the relationship between the exceedance probability of streamflow and its magnitude. This curve is used widely for flood risk analysis, water quality management, and the design of hydroelectric power plants (among others). Several mathematical formulations have been proposed to mimic the FDC. Yet, these functions are often not flexible enough to portray accurately the functional shape of the FDC for a large range of catchments. \cite{vrugt2013} introduced the soil water characteristic (SWC) of van Genuchten \citep{genuchten1980} as new parametric expression of the FDC for diagnostic model evaluation with DREAM$_{\rm (ABC)}$. \cite{sadegh2016} build on the work of \cite{vrugt2013} and compared several models of the SWC against their counterparts published in the literature. These new expressions were shown to fit well the empirical FDCs of the 438 watersheds of the MOPEX data set. Here, we present a MATLAB toolbox, called FDCFIT which contains the fifteen different FDC functions described in \cite{sadegh2016} and returns the values of their coefficients for a given discharge record, along with graphical output of the fit. Two case studies are used to illustrate the main capabilities and functionalities of the FDCFIT toolbox.

## Getting Started

### Installing

* Download and unzip the zip file 'MATLAB_pcode_FDCFIT_V1.0.zip'.
* Add the toolbox to your MATLAB search path by running the script 'install_FDCFIT.m' available in 'MATLAB_pcode_FDCFIT_V1.0'
* You are ready to run the examples.

### Executing program

* After intalling, you can simply direct to each example folder and execute the local 'example_X.m' file.
* Please make sure you read carefully the instructions (i.e., green comments) in 'install_MODELAVG.m' and the manual in PDF !!!  
* Please also refer to the attached paper by Sadegh, Vrugt, Gupta and Xu (2016) 

## Authors

* Vrugt, Jasper A. (jasper@uci.edu) 

## Version History

* 1.0
    * Initial Release

## Acknowledgments
