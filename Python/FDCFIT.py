# ####################################################################### #
#                                                                         #
#  FFFFFFFFFFF DDDDDDDDD  CCCCCCCCCC FFFFFFFFFFF IIIIIIIIIIII TTTTTTTTTT  #
#  FFFFFFFFFFF DDDDDDDDDD CCCCCCCCC  FFFFFFFFFFF  IIIIIIIIII  TTTTTTTTTT  #
#  FF          DD      DD CC         FF               II          TT      #
#  FF          DD      DD CC         FF               II          TT      #
#  FFFFFF      DD      DD CC         FFFFFF           II          TT      #
#  FF          DD      DD CC         FF               II          TT      #
#  FF          DDDDDDDDDD CCCCCCCCC  FF           IIIIIIIIII      TT      #
#  FF          DDDDDDDDD  CCCCCCCCCC FF          IIIIIIIIIIII     TT      #
#                                                                         #
# ####################################################################### #
#                                                                         #
# This code determines the fitting coefficients of a large suite of       #
# parametric expressions of the flow duration curve. This includes a new  #
# class of functions based on the water retention functions of van        #
# Genuchten (1980) and Kosugi (1996).                                     #
#                                                                         #
###########################################################################
#                                                                         #
# SYNOPSIS: [x,RMSE] = FDCFIT(FDCPar,model,e,y,method)                    #
#           [x,RMSE] = FDCFIT(FDCPar,model,e,y,method,options)            #
#  where                                                                  #
#   FDCPar    [input] Structure with settings for FDCFIT                  #
#    .d           # unknown parameters of parametric FDC expression       #
#    .form        Forward or inverse formulation of the FDC               #
#      = 'e'      Forward expression: e_s = f(y)                          #
#      = 'y'      Backward/inverse expression: y_s = f(e)                 #
#    .n           # discharge measurements [= numel(y)]  return by FDCFIT #
#    .model       Full name of FDC expression            return by FDCFIT #
#    .acronym     Acroynym FDC expression                return by FDCFIT #
#   model     [input] Integer of parametric FDC expression                #
#      = 1          Full name: lognormal-2         Acronym: LN-2          #
#      = 2          Full name: gumbel              Acronym: G             #
#      = 3          Full name: logistic            Acronym: LG            #
#      = 4          Full name: logarithmic         Acronym: LOG           #
#      = 5          Full name: power               Acronym: PW            #
#      = 6          Full name: quimpo              Acronym: Q             #
#      = 7          Full name: viola               Acronym: V             #
#      = 8          Full name: genuchten-2         Acronym: VG-2          #
#      = 9          Full name: kosugi-2            Acronym: K-2           #
#      = 10         Full name: lognormal-3         Acronym: LN-3          #
#      = 11         Full name: pareto              Acronym: GP            #
#      = 12         Full name: gev                 Acronym: GEV           #
#      = 13         Full name: franchini           Acronym: F             #
#      = 14         Full name: genuchten-3         Acronym: VG-3          #
#      = 15         Full name: kosugi-3            Acronym: K-3           #
#   e         [input] nx1 vector of empirical exceedance probabilities    #
#   y         [input] nx1 vector of empirical (= measured) discharges     #
#   method    [input] Name (string) of parameter estimation method        #
#    = 'lm'         Levenberg-Marquardt method                            #
#    = 'sp'         Nelder-Mead simplex algorithm                         #
#    = 'de'         Differential evolution                                #
#    = 'cma'        Covariance-matrix adaptation                          #
#   options   [input] Optional structure for algorithmic variables        #
#    .N             # trials 'lm' and 'sp' methods        DEF: 5          #
#    .P             Population size for 'de' and 'cma'    DEF: 25         #
#    .CR            Crossover value for 'de' algorithm    DEF: 0.8        #
#    .TolX          Error tolerance on parameter values   DEF: 1e-3       #
#    .TolFun        Error tolerance on objective function DEF: 1e-4       #
#    .MaxFunEvals   Maximum # function evaluations        DEF: 1e4        #
#    .type          Type of FDC derived from daily discharge data, y      #
#     = 'day'       Daily flow duration curve                             #
#     = 'week'      Weekly flow duration curve                            #
#     = 'month'     Monthly flow duration curve                           #
#     = 'annual'    Annual flow duration curve                            #
#    .print         Output writing screen (tables/figs)   DEF: 'yes'      #
#     = 'no'        No output writing to the screen                       #
#     = 'yes'       Output writing to the screen          DEFault         #
#   x         [outpt] Maximum likelihood parameters FDC expressions       #
#   RMSE      [outpt] Root Mean Square Error of FDC expression            #
#                                                                         #
# ####################################################################### #
#                                                                         #
#  BUILT-IN CASE STUDIES:                                                 #
#   Example 1   French Broad River - wettest of MOPEX data set            #
#   Example 2   Guadalupe River basin - driest of MOPEX data set          #
#   Example 3   Another 10 watersheds of the MOPEX data set               #
#                                                                         #
# ####################################################################### #
#                                                                         #
# LITERATURE:                                                             #
#  Vrugt, J.A. (2017), FDCFIT: A MATLAB toolbox of parametric             #
#      expressions of the Flow duration curve, Manual, pp. 1-37.          #
#  Sadegh, M., J.A. Vrugt, X. Cu, and H.V. Gupta, (2016), The soil water  #
#      characteristic as new class of parametric expressions of the flow  #
#      duration curve, Journal of Hydrology, 535, pp. 438-456,            #
#          https://doi.org/10.1016/j.jhydrol.2016.01.027                  #
#  Vrugt, J.A., and M. Sadegh (2013), Toward diagnostic model             #
#      calibration and evaluation: Approximate Bayesian computation,      #
#      Water Resources Research, 49, 4335–4345,                           #
#          https://doi.org/10.1002/wrcr.20354                             #
#                                                                         #
# ####################################################################### #
#                                                                         #
#  PYTHON CODE:                                                           #
#  © Written by Jasper A. Vrugt using GPT-4 OpenAI's language model       # 
#    University of California Irvine                                      #
#  Version 2.0    Dec 2024                                                #
#                                                                         #
# ####################################################################### #

import numpy as np                                      
import os, sys
from scipy.optimize import least_squares, minimize 

example_dir = os.getcwd()					                    # Add current directory to Python path
if example_dir not in sys.path:
    sys.path.append(example_dir)

parent_dir = os.path.abspath(os.path.join(example_dir, '..'))   # Go up one directory
sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	    # Add miscellaneous directory to Python path
from FDCFIT_functions import *				                    # Import functions

def FDCFIT(FDCPar, model, e, y, method, options=None):
    # ####################################################################### #    
    # This function determines the fitting coefficients of a large suite of   #
    # models for the flow duration curve (FDC).                               #
    # ####################################################################### #    

    if options is None:
        options = {}

    # Check if method is valid and make sure model is integer, JAV Python
    method = method.lower(); model = int(model)

    # Check setup of user inputs and validate
    FDCPar, options, f_names = FDCFIT_check(FDCPar, model, e, y, method, options)

    # Setup problem (parameter ranges, model, etc.)
    X, FDCPar, Par_info, options, lsq_options, n_crash, units, math_str = FDCFIT_setup(FDCPar, model, e, method, options, f_names)

    # Dynamic part: Optimize based on method
    if method == 'lm':  # Levenberg-Marquardt method
        err = np.full(options['N'], np.nan); it_MAP = np.full(options['N'], np.nan) 
        MAP = np.full((options['N'], FDCPar['d']), np.nan)
        for j in range(options['N']):
            try:
                result = least_squares(FDCFIT_functions, X[j, :FDCPar['d']], args=(FDCPar, e, y, 2), bounds=(Par_info['min'], Par_info['max']), method='dogbox') #, **lsq_options)
                MAP[j, :] = result.x
                err[j] = result.cost
                it_MAP[j] = result.nfev
            except Exception:
                n_crash += 1
        RMSE_MAP = np.sqrt(err / FDCPar['n'])

    elif method == 'sp':  # Nelder-Mead Simplex method
        err = np.full(options['N'], np.nan); it_MAP = np.full(options['N'], np.nan) 
        MAP = np.full((options['N'], FDCPar['d']), np.nan)
        for j in range(options['N']):
            try:
                res = minimize(FDCFIT_functions, X[j, :FDCPar['d']], args=(FDCPar, e, y, 1), method='Nelder-Mead', options=lsq_options)
                MAP[j, :] = res.x
                err[j] = res.fun
                it_MAP[j] = res.nfev
            except Exception:
                n_crash += 1
        RMSE_MAP = np.sqrt(err / FDCPar['n'])

    elif method == 'de':  # Differential evolution
        MAP, RMSE_MAP, it_MAP = de_code(X, FDCPar, e, y, Par_info, options)

    elif method == 'cma':  # Covariance-matrix adaptation
        MAP, RMSE_MAP, it_MAP = cmaes_code(FDCPar, e, y, options)

    else:
        raise ValueError("FDCFIT: Unknown parameter estimation method")

    print('\n')

    # If there were crashes, log them
    if n_crash > 0:
        with open('warning_file.txt', 'a+') as fid:
            evalstr = f"FDCFIT WARNING: Of the {options['N']} trials, {n_crash} were unsuccessful (crashed)\n"
            fid.write(evalstr)
            print(evalstr)
            
            evalstr = ('FDCFIT SUGGESTION: If too many LM trials crashed: Use global search '
                       'with Differential Evolution (= DE) or CMAES (= CMA) instead; or try '
                       'Nelder-Mead Simplex (= SP)\n')
            fid.write(evalstr)
            print(evalstr)

    # Final processing and output
    map, RMSE_map, it_map, str_ = FDCFIT_end(FDCPar, MAP, RMSE_MAP, it_MAP)

    if options.get('print', 'yes') == 'yes':
        FDCFIT_plot(FDCPar, e, y, method, map, RMSE_map, str_, units, math_str)

    return map, RMSE_map, it_map