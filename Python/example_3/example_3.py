###########################################################################
#  FFFFFFFFFFF DDDDDDDDD  CCCCCCCCCC FFFFFFFFFFF IIIIIIIIIIII TTTTTTTTTT  #
#  FFFFFFFFFFF DDDDDDDDDD CCCCCCCCC  FFFFFFFFFFF  IIIIIIIIII  TTTTTTTTTT  #
#  FF          DD      DD CC         FF               II          TT      #
#  FF          DD      DD CC         FF               II          TT      #
#  FFFFFF      DD      DD CC         FFFFFF           II          TT      #
#  FF          DD      DD CC         FF               II          TT      #
#  FF          DDDDDDDDDD CCCCCCCCC  FF           IIIIIIIIII      TT      #
#  FF          DDDDDDDDD  CCCCCCCCCC FF          IIIIIIIIIIII     TT      #
###########################################################################

import numpy as np
import os, sys

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from FDCFIT import FDCFIT

sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	    # Add miscellaneous directory to Python path
from calc_FDC import calc_FDC	                                # Import functions

# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      333333   #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE              33   #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE          333   #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE              33   #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE      333333   #
#                                                                         #
# ####################################################################### #

# Following models are part of the FDCFIT toolbox
# model       name      Acronym used in paper Sadegh, Vrugt et al., 2015
#   1   'lognormal-2'   LN-2
#   2   'gumbel'        G
#   3   'logistic'      LG
#   4   'logarithmic'   LOG
#   5   'power'         PW
#   6   'quimpo'        Q
#   7   'viola'         V
#   8   'genuchten-2'   VG-2
#   9   'kosugi-2'      K-2
#   10  'lognormal-3'   LN-3
#   11  'pareto'        GP
#   12  'gev'           GEV
#   13  'franchini'     F
#   14  'genuchten-3'   VG-3
#   15  'kosugi-3'      K-3

# Define parametric expression for FDC
model = 14
# Define the FDCPar structure (parameters for the FDC fit)
FDCPar = {'form': 'e'}
method = 'LM';      # Optimization method 
                    # 'LM'    Levenberg Marquardt          (= local method)
                    # 'SP'    Nelder-Mead Simplex          (= local method)
                    # 'CMA'   Covariance Matrix Adaptation (= global methd)
                    # 'DE'    Differential Evolution       (= global methd)

# Define the options dictionary
options = {
    'N': 5,                 # Number of trials with optimization method
    'TolX': 1e-2,           # Termination criteria on parameters
    'TolFun': 1e-3,         # Termination criteria objective function (SSE)
    'MaxFunEvals': 1e4,     # Maximum # function evaluations each trial
    'type': 'day',          # Time scale of flow duration curve
    'print': 'yes'          # Output to screen ('yes' or 'no')
}

# Example 3: Load the discharge data of basin j of MOPEX data set
# (unzip file "MOPEX_remaining_watersheds.rar" for all MOPEX data sets)
watershed_no = 10

file_name = f"data_{watershed_no}.txt"              # Construct the file name dynamically
data = np.loadtxt(file_name)                        # Load the data from the text file
ii = (data == 0)				                    # Zero values of discharge
p_0 = np.sum(ii) / data.size  		                # Probability of zero flows
Y = data[(data>0)]      			                # Remaining data excluding zeros
Y = Y[~np.isnan(Y)] 				                # Remove nan from discharge series
N = Y.size					                        # How many values of Y?
y_s = np.sort(Y)				                    # Sort the discharge data in ascending order
e_n = 1 - ((np.arange(1, N+1)) - 0.5) / N	        # Calculate the exceedance probabilities
# e, y, p_0 = calc_FDC(data, options['type'], 1)    # Call FDC function for dly files [MOPEX]

if __name__ == '__main__':
	# Run FDCFIT (fit the curve using the chosen model)
	map, RMSE_map, it_map = FDCFIT(FDCPar, model, e_n, y_s, method, options)
