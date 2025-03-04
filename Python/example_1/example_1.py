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
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE        1111   #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE           11 11   #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       11  11   #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE              11   #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE          11   #
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

# Define the model selection (choose the model between 1 and 15)
model = 10

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

# Example 1: French Broad River - wettest of MOPEX data set
ID_watershed = '08167500'

# Call FDC function
e, y, p_0 = calc_FDC(ID_watershed, options['type'], 6)

if __name__ == '__main__':
	# Run FDCFIT (fit the curve using the chosen model)
	map, RMSE_map, it_map = FDCFIT(FDCPar, model, e, y, method, options)

