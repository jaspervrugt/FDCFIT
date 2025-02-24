# FDCFIT: Fitting of flow duration curves using parametric formulations in MATLAB and Python

## Description

The flow duration curve (FDC) is a signature catchment characteristic that depicts graphically the relationship between the exceedance probability of streamflow and its magnitude. This curve is used widely for flood risk analysis, water quality management, and the design of hydroelectric power plants (among others). Several mathematical formulations have been proposed to mimic the FDC. Yet, these functions are often not flexible enough to portray accurately the functional shape of the FDC for a large range of catchments. _Vrugt and Sadegh_ (2013) introduced the soil water characteristic (SWC) of _van Genuchten_ (1980) as new parametric expression of the FDC for diagnostic model evaluation with the DREAM$_{(ABC)}$ algorithm. _Sadegh et al._ (2016) build on the work of _Vrugt and Sadegh_ (2013) and compared several models of the SWC against their counterparts published in the literature. These new expressions were shown to fit well the empirical FDCs of the 438 watersheds of the MOPEX data set. Here, we present a MATLAB and Python toolbox, called FDCFIT which implements 19 different two- and three-parameter expressions for the FDC (15 from _Sadegh et al_., 2016) and returns optimized values of their coefficients for the discharge record at hand, along with graphical output of measured and fitted FDC's in linear and logarithmic space. Users can choice between local (Levenberg-Marquardt & Nelder-Mead Simplex) and population-based (Covariance Matrix Adaption & Differential Evolution) search algorithms. Three built-in case studies illustrate the main capabilities and functionalities of the FDCFIT toolbox.

## Getting Started

### Installing: MATLAB

* Download and unzip the zip file 'MATLAB_code_FDCFIT_V2.0.zip' in a directory 'FDCFIT'
* Add the toolbox to your MATLAB search path by running the script 'install_FDCFIT.m' available in the root directory
* You are ready to run the examples

### Executing program

* After intalling, you can simply direct to each example folder and execute the local 'example_X.m' file
* Please make sure you read carefully the instructions (i.e., green comments) in 'install_FDCFIT.m' and the manual !!!  

### Installing: Python

* Download and unzip the zip file 'Python_code_FDCFIT_V2.0.zip' to a directory called 'FDCFIT'

### Executing program

* Go to Command Prompt and directory of example_X in the root of FDCFIT
* Now you can execute this example by typing "python example_X.py".
* Instructions can be found in the file 'FDCFIT.py' and in the manual !!!
  
## Authors

* Vrugt, Jasper A. (jasper@uci.edu) 

## Literature

1. Vrugt, J.A. (2017), FDCFIT: A MATLAB toolbox of parametric expressions for the flow duration curve, Manual, pp. 1-37
2. Sadegh, M., J.A. Vrugt, X. Cu, and H.V. Gupta, (2016), The soil water characteristic as new class of parametric expressions of the flow duration curve, _Journal of Hydrology_, 535, pp. 438-456, https://doi.org/10.1016/j.jhydrol.2016.01.027
3. Vrugt, J.A., and M. Sadegh (2013), Toward diagnostic model calibration and evaluation: Approximate Bayesian computation, _Water Resources Research_, 49, 4335–4345, https://doi.org/10.1002/wrcr.20354

## Version History

* 1.0
    * Initial Release
* 2.0
    * Cleaning up code and enhancing postprocessing capabilities
    * Python implementation

## Built-in Case Studies

1. Example 1: French Broad River - wettest of MOPEX data set
2. Example 2: Guadalupe River basin - driest of MOPEX data set
3. Example 3: Ten different MOPEX data sets

## Acknowledgments
