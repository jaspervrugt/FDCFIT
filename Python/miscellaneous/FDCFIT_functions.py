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
import os, sys, math
from scipy.optimize import least_squares, minimize
from scipy.stats import gamma
from scipy.special import erfc, erfcinv
import scipy.stats as stats
import logging
import matplotlib.pyplot as plt
from matplotlib import cm
from pandas.plotting import table
from matplotlib.backends.backend_pdf import PdfPages
from screeninfo import get_monitors


def FDCFIT_check(FDCPar, model, e, y, method, options):
    # ####################################################################### #
    #                                                                         #
    # Checks the settings defined by the user                                 #
    #                                                                         #
    # SYNOPSIS: [FDCPar,options,f_names] = FDCFIT_check(FDCPar,model,e,y, ... #
    #               method,options)                                           #
    #  where                                                                  #
    #   FDCPar    [input] Structure with settings for FDCFIT                  #
    #   model     [input] Integer of parametric FDC expression                #
    #      = 1          lognormal-2'                                          #
    #      = 2          gumbel                                                #
    #      = 3          logistic                                              #
    #      = 4          logarithmic                                           #
    #      = 5          power                                                 #
    #      = 6          quimpo                                                #
    #      = 7          viola                                                 #
    #      = 8          genuchten-2                                           #
    #      = 9          kosugi-2                                              #
    #      = 10         lognormal-3                                           #
    #      = 11         pareto                                                #
    #      = 12         gev                                                   #
    #      = 13         franchini                                             #
    #      = 14         genuchten-3                                           #
    #      = 15         kosugi-3                                              #
    #      = 16         weibull                                               #
    #      = 17         non_central_F                                         #
    #      = 18         exponential                                           #
    #      = 19         gamma                                                 #
    #   e         [input] nx1 vector of empirical exceedance probabilities    #
    #   y         [input] nx1 vector of empirical (= measured) discharges     #
    #   method    [input] Name (string) of parameter estimation method        #
    #    = 'lm'         Levenberg-Marquardt method                            #
    #    = 'sp'         Nelder-Mead simplex algorithm                         #
    #    = 'de'         Differential evolution                                #
    #    = 'cma'        Covariance-matrix adaptation                          #
    #   options   [input] Optional structure for algorithmic variables        #
    #   FDCPar    [outpt] Structure with settings for FDCFIT [verified]       #
    #   options   [outpt] Optional structure algorithmic variables [verified] #
    #   f_names   [outpt] Name (string) of suite FDC parametric expressions   #
    #                                                                         #
    # (c) Written by Jasper A. Vrugt, Dec. 2014                               #
    # University of California Irvine                                         #
    #                                                                         #
    # ####################################################################### #
        
    # Define model names
    f_names = [
        'lognormal-2', 'gumbel', 'logistic', 'logarithmic', 'power', 'quimpo',
        'viola', 'genuchten-2', 'kosugi-2', 'lognormal-3', 'pareto', 'gev',
        'franchini', 'genuchten-3', 'kosugi-3', 'weibull', 'non_central_F', 'exponential', 'gamma'
    ]
    logger = logging.getLogger(__name__)  # Creates a logger with the name of the current module

    # Print header information
    print('  ------------------------------------------------------------------             ')
    print('  FFFFFFFFF  DDDDDDDD   CCCCCCCCC  FFFFFFFFF  IIIIIIIII  TTTTTTTTTTT             ')
    print('  FFFFFFFF   DDDDDDDDD  CCCCCCCC   FFFFFFFF   IIIIIIIII  TTTTTTTTTTT             ')
    print('  FFF        DDD   DDD  CCC        FFF           III     TT  TTT  TT             ')
    print('  FFF        DDD   DDD  CCC        FFF           III     TT  TTT  TT             ')
    print('  FFFFFF     DDD   DDD  CCC        FFFFFF        III         TTT         /^ ^\   ')
    print('  FFFFFF     DDD   DDD  CCC        FFFFFF        III         TTT        / 0 0 \  ')
    print('  FFF        DDD   DDD  CCC        FFF           III         TTT        V\ Y /V  ')
    print('  FFF        DDD   DDD  CCC        FFF           III         TTT         / - \   ')
    print('  FFF        DDDDDDDDD  CCCCCCCC   FFF       IIIIIIIIII      TTT        /     |  ')
    print('  FFF        DDDDDDDD   CCCCCCCCC  FFF       IIIIIIIIII      TTT        V__) ||  ')
    print('  ------------------------------------------------------------------ ')
    print('  © Jasper A. Vrugt, University of California Irvine & GPT-4 OpenAI''s language model')
    print('    ________________________________________________________________________')
    print('    Version 2.0, Dec. 2024, Beta-release: MATLAB implementation is benchmark')
    print('\n')

    # Open warning file for output
    with open('warning_file.txt', 'w+', encoding='utf-8') as fid:
        fid.write('-------------- FDCFIT warning file --------------\n')

        # Check FDCPar input argument
        if FDCPar is None or not isinstance(FDCPar, dict):
            raise ValueError('FDCFIT ERROR: FDCPar should be a non-empty dictionary')

        # Check model input argument
        if model is None or not isinstance(model, int) or not (1 <= model <= 19):
            raise ValueError(f'FDCFIT ERROR: model should be an integer between 1 and 19, got {model}')
 
        # Check e and y
        if e is None or not isinstance(e, np.ndarray) or e.ndim != 1:
            raise ValueError('FDCFIT ERROR: e should be a 1D array (nx1 vector) of exceedance probabilities')

        if y is None or not isinstance(y, np.ndarray) or y.ndim != 1:
            raise ValueError('FDCFIT ERROR: y should be a 1D array (nx1 vector) of discharges')

        if len(e) != len(y):
            raise ValueError('FDCFIT ERROR: e and y must have the same length')

        if len(e) < 3:
            raise ValueError('FDCFIT ERROR: Insufficient number of observations for e and y')

        if len(e) < 25:
            warning_msg = 'FDCFIT WARNING: Rather small number of observations of discharge/exceedance probability\n'
            print(warning_msg)
            fid.write(warning_msg)

        # Check options input argument
        if options is None or not isinstance(options, dict):
            raise ValueError('FDCFIT ERROR: options should be a non-empty dictionary')

        # Check FDCPar.form field
        if 'form' not in FDCPar:
            raise ValueError('FDCFIT ERROR: The field "form" in FDCPar is not defined. Set it to "e" or "y"')
        
        if not isinstance(FDCPar['form'], str):
            raise ValueError('FDCFIT ERROR: "form" should be a string ("e" or "y")')

        FDCPar['form'] = FDCPar['form'].lower()
        if FDCPar['form'] not in ['e', 'y']:
            raise ValueError('FDCFIT ERROR: "form" should be "e" or "y"')

        # Check method input argument
        if not isinstance(method, str) or method.lower() not in ['lm', 'sp', 'de', 'cma']:
            raise ValueError("FDCFIT ERROR: Unknown optimization method. Use 'LM', 'SP', 'DE' or 'CMA'")

        # Check options.N field for LM/SP methods
        if method.lower() in ['lm', 'sp']:
            if 'N' not in options:
                warning_msg = "FDCFIT WARNING: 'N' (number of trials with SP/LM) is not defined. Using default N=5\n"
                print(warning_msg)
                fid.write(warning_msg)
                options['N'] = int(5)  # Set default value
                
            elif not isinstance(options['N'], int): 
                raise ValueError("FDCFIT ERROR: 'N' should be an integer larger than zero")
            
            options['N'] = int(options['N'])
            if options['N'] <= 0:
                raise ValueError("FDCFIT ERROR: 'N' should be a positive integer")
            if options['N'] < 5:
                warning_msg = f'FDCFIT WARNING: The value of N ({options["N"]}) is rather small. Default is N=5\n'
                print(warning_msg)
                fid.write(warning_msg)
            elif options['N'] > 25:
                warning_msg = f'FDCFIT WARNING: The value of N ({options["N"]}) is quite large. Default is N=5\n'
                print(warning_msg)
                fid.write(warning_msg)

            if options['N'] % 1 != 0:
                warning_msg = f"FDCFIT WARNING: 'N' should be an integer. Using default N=5\n"
                print(warning_msg)
                fid.write(warning_msg)
                options['N'] = int(5)
            options['N'] = int(options['N'])
        # Default values
        if 'TolX' not in options:
            logger.warning("FDCFIT WARNING: Field 'TolX' of structure options not defined by user --> options.TolX = 1e-2 (default value)")
            options['TolX'] = 1e-2
        else:
            if options['TolX'] is None:
                logger.warning("FDCFIT WARNING: Field 'TolX' of structure options is left empty by user --> default settings options.TolX = 1e-2")
                options['TolX'] = 1e-2
            elif not isinstance(options['TolX'], (int, float)):
                logger.warning("FDCFIT WARNING: The field 'TolX' of structure options should be a numerical value --> resort to default setting options.TolX = 1e-2")
                options['TolX'] = 1e-2
            elif options['TolX'] < 0:
                logger.warning("FDCFIT WARNING: The field 'TolX' of structure options cannot be negative (default is options.TolX = 1e-2)")
                options['TolX'] = 1e-2
            elif options['TolX'] > 0.1:
                logger.warning(f"FDCFIT WARNING: The value of {options['TolX']} for options.TolX is rather large (default is options.TolX = 1e-2)")

        # Check TolFun
        if 'TolFun' not in options:
            logger.warning("FDCFIT WARNING: Field 'TolFun' of structure options not defined by user --> options.TolFun = 1e-4 (default value)")
            options['TolFun'] = 1e-4
        else:
            if options['TolFun'] is None:
                logger.warning("FDCFIT WARNING: Field 'TolFun' of structure options is left empty by user --> default settings options.TolFun = 1e-4")
                options['TolFun'] = 1e-4
            elif not isinstance(options['TolFun'], (int, float)):
                logger.warning("FDCFIT WARNING: The field 'TolFun' of structure options should be a numerical value --> resort to default setting options.TolFun = 1e-4")
                options['TolFun'] = 1e-4
            elif options['TolFun'] < 0:
                logger.warning("FDCFIT WARNING: The field 'TolFun' of structure options cannot be negative (default is options.TolFun = 1e-4)")
                options['TolFun'] = 1e-4
            elif options['TolFun'] > 0.1:
                logger.warning(f"FDCFIT WARNING: The value of {options['TolFun']} for options.TolFun is rather large (default is options.TolFun = 1e-4)")

        # Check MaxFunEvals
        if 'MaxFunEvals' not in options:
            logger.warning("FDCFIT WARNING: Field 'MaxFunEvals' of structure options not defined by user --> options.MaxFunEvals = 1e4 (default value)")
            options['MaxFunEvals'] = int(1e4)
        else:
            if options['MaxFunEvals'] is None:
                logger.warning("FDCFIT WARNING: Field 'MaxFunEvals' of structure options is left empty by user --> default settings options.MaxFunEvals = 1e4")
                options['MaxFunEvals'] = int(1e4)
            elif not isinstance(options['MaxFunEvals'], (int, float)):
                logger.warning("FDCFIT WARNING: The field 'MaxFunEvals' of structure options should be a numerical value --> resort to default setting options.MaxFunEvals = 1e4")
                options['MaxFunEvals'] = int(1e4)
            elif options['MaxFunEvals'] < 0:
                logger.warning("FDCFIT WARNING: The field 'MaxFunEvals' of structure options cannot be negative (default is options.MaxFunEvals = 1e4)")
                options['MaxFunEvals'] = int(1e4)
            elif options['MaxFunEvals'] < 51:
                raise ValueError("FDCFIT ERROR: The field 'MaxFunEvals' of structure options is set too small to provide reliable results (default is options.MaxFunEvals = 1e4)")
            elif options['MaxFunEvals'] > 1e5:
                logger.warning(f"FDCFIT WARNING: The value of {options['MaxFunEvals']} for optim.MaxFunEvals is rather large (default is options.MaxFunEvals = 1e4)")
            options['MaxFunEvals'] = int(options['MaxFunEvals'])
        # Check 'type'
        if 'type' in options:
            if options['type'] is None:
                logger.warning("FDCFIT WARNING: The field 'type' of structure options should not be empty but list a string (one of 'day'/'week'/'month'/'annual')")
            elif not isinstance(options['type'], str):
                logger.warning("FDCFIT WARNING: The field 'type' of structure options should be a string (one of 'day'/'week'/'month'/'annual')")
            elif options['type'] not in ['day', 'week', 'month', 'annual']:
                logger.warning("FDCFIT WARNING: The field 'type' of structure options should be a string (one of 'day'/'week'/'month'/'annual')")

        # Check 'print'
        if 'print' in options:
            if options['print'] is None:
                logger.warning("FDCFIT WARNING: The field 'print' of structure options should not be empty but list a string with either 'yes' or 'no' (resort to default 'yes')")
                options['print'] = 'yes'
            elif not isinstance(options['print'], str):
                logger.warning("FDCFIT WARNING: The field 'print' of structure options should be a string with either 'yes' or 'no' (resort to default 'yes')")
                options['print'] = 'yes'
            elif options['print'] not in ['yes', 'no']:
                logger.warning("FDCFIT WARNING: The field 'print' of structure options should be a string with either 'yes' or 'no' (resort to default 'yes')")
                options['print'] = 'yes'

        # DE/CMAES-specific checks
        if method.lower() in ['de', 'cma']:
            if 'P' not in options:
                logger.warning("FDCFIT WARNING: Field 'P' (population size) of structure options not defined by user --> options.P = 25 (default value)")
                options['P'] = int(25)
            elif options['P'] is None:
                logger.warning("FDCFIT WARNING: Field 'P' of structure options is left empty by user --> default settings options.P = 25")
                options['P'] = int(25)
            elif not isinstance(options['P'], (int, float)):
                logger.warning("FDCFIT WARNING: The field 'P' of structure options should be a numerical value --> resort to default setting options.P = 25")
                options['P'] = int(25)
            elif options['P'] < 0:
                logger.warning("FDCFIT WARNING: The field 'P' of structure options cannot be negative: Set to at least 'P = 25' individuals for global search with " + method.upper())
                options['P'] = int(25)
            elif options['P'] < 10:
                logger.warning(f"FDCFIT WARNING: The field 'P' of structure options is set rather small: At least 'P = 25' individuals for global search with {method.upper()}")
            elif options['P'] > 50:
                logger.warning(f"FDCFIT WARNING: The field 'P' of structure options is set rather large: A population of 'P = 50' should be sufficient for estimating a few parameters with {method.upper()}")
            options['P'] = int(options['P'])
 
        # DE-specific checks
        if method.lower() == 'de':
            if 'CR' not in options:
                logger.warning("FDCFIT WARNING: Field 'CR' (crossover value) of structure DE not defined by user --> options.CR = 0.9 (default value)")
                options['CR'] = 0.9
            elif options['CR'] is None:
                logger.warning("FDCFIT WARNING: Field 'CR' of structure options is left empty by user --> default settings options.CR = 0.9")
                options['CR'] = 0.9
            elif not isinstance(options['CR'], (int, float)):
                logger.warning("FDCFIT WARNING: The field 'CR' of structure options should be a numerical value --> resort to default setting options.CR = 0.9")
                options['CR'] = 0.9
            elif options['CR'] < 0:
                logger.warning("FDCFIT WARNING: The field 'CR' of structure options cannot be smaller than zero (default is options.CR = 0.9)")
                options['CR'] = 0.9
            elif options['CR'] > 1:
                logger.warning("FDCFIT ERROR: The field 'CR' of structure options cannot be larger than one (default is options.CR = 0.9)")
                options['CR'] = 0.9
            elif options['CR'] < 0.4:
                logger.warning("FDCFIT WARNING: The field 'CR' of structure options cannot be smaller than zero (default is options.CR = 0.6)")
                options['CR'] = 0.6
            elif options['CR'] > 0.9:
                logger.warning("FDCFIT ERROR: The field 'CR' of structure options is set rather large (default is options.CR = 0.9)")
                options['CR'] = 0.9

    return FDCPar, options, f_names


def FDCFIT_setup(FDCPar, model, e, method, options, f_names):
    # ####################################################################### #
    #                                                                         #
    # Defines parameter ranges for each model class and formulation           #
    #                                                                         #
    # SYNOPSIS: [X,FDCPar,Par_info,options,lsq_options,index,units, ...       #
    #               math_str,count] = FDCFIT_setup(FDCPar,model,e,method, ... #
    #               options,f_names)                                          #
    #  where                                                                  #
    #   FDCPar    [input] Structure with settings for FDCFIT                  #
    #   model     [input] Integer of parametric FDC expression                #
    #   e         [input] nx1 vector of empirical exceedance probabilities    #
    #   method    [input] Name (string) of parameter estimation method        #
    #   options   [input] Optional structure for algorithmic variables        #
    #   f_names   [input] Name (string) of suite FDC parametric expressions   #
    #   X         [outpt] initial draws of parameter values                   #
    #   FDCPar    [outpt] Dictionary with algorithmic parameters              #
    #   Par_info  [outpt] Dictionary with parameter information               #
    #   options   [outpt] Optional structure algorithmic variables [verified] #
    #   lsq_options [otp] Dictionary with options least squares method        # 
    #   n_crash   [outpt] Initialize at zero number of trials that crashed    #   
    #   units     [outpt] Unit specification                                  # 
    #   math_str  [outpt] String with mathematical function model formulation #
    #                                                                         #
    # (c) Written by Jasper A. Vrugt, Dec. 2014                               #
    # University of California Irvine                                         #
    #                                                                         #
    # ####################################################################### #

    # Initialize variables
    n_crash = 0
    logging.basicConfig(filename = 'warning_file.txt', level = logging.WARNING)

    # conversion to Python, JAV
    model = int(model) - 1

    # Define model acronyms
    model_acronyms = ['LN', 'G', 'LG', 'LOG', 'PW', 'Q', 'V', 'VG', 'K', 'LN', 'GP', 'GEV', 'FS', 'VG', 'K', 'WBL', 'NCF', 'EXP', 'GAM']
    
    # Define the number of parameters for each model - match acronyms
    n_pars = [2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 2, 3, 2, 2]
    
    # Determine the number of parameters for the selected model
    FDCPar['d'] = n_pars[model]
    
    # Store model details
    FDCPar['acronym'] = model_acronyms[model]
    FDCPar['model'] = f_names[model]
    FDCPar['n'] = len(e)

    # Define mathematical expressions and parameter ranges for each model
    if FDCPar['model'] == 'lognormal-2':
        min_vals = np.array([-np.inf, 0])
        max_vals = np.array([np.inf, np.inf])
        units = ['(L/T)', '(L/T)']
        if FDCPar['form'] == 'y':
            math_str = r'${\rm Lognormal-2:} y = \exp\left(a_{\rm LN} - \sqrt{2}b_{\rm LN} \text{erfc}^{-1}\bigl(2(1 - e_{\rm n})\bigr)\right)$'
        else:
            math_str = r'${\rm Lognormal-2:} e_{\rm n} = 1-\frac{1}{2}\text{erfc} \left(\frac{a_{\rm LN} - \log(y)}{\sqrt{2} \: b_{\rm LN}}\right)$'
    
    elif FDCPar['model'] == 'gumbel':
        min_vals = np.array([-np.inf, 0])
        max_vals = np.array([np.inf, np.inf])
        units = ['(L/T)', '(L/T)']
        if FDCPar['form'] == 'y':
            math_str = r'${\rm Gumbel:} y = a_{\rm G} - b_{\rm G}\log\left[\log\left(\frac{1}{(1 - e_{\rm n})}\right) \right]$'
        else:
            math_str = r'${\rm Gumbel:} e_{\rm n} = 1 - \exp\left[-\exp\left(\frac{a_{\rm G} - y}{b_{\rm G}}\right)\right]$'
    
    elif FDCPar['model'] == 'logistic':
        min_vals = np.array([-np.inf, 0])
        max_vals = np.array([np.inf, np.inf])
        units = ['(L/T)', '(L/T)']
        if FDCPar['form'] == 'y':
            math_str = r'${\rm Logistic:} y = a_{\rm LG} - b_{\rm LG}\log\left(\frac{1}{(1-e_{\rm n})} - 1\right)$'
        else:
            math_str = r'${\rm Logistic:} e_{\rm n} = 1 - \left[1 + \exp\left(\frac{a_{\rm LG} - y}{b_{\rm LG}}\right)\right]^{-1}$'

    elif FDCPar['model'] == 'logarithmic':
        min_vals = np.array([-np.inf, -np.inf])
        max_vals = np.array([0, np.inf])
        units = ['(L/T)', '(L/T)']
        if FDCPar['form'] == 'y':
            math_str = r'${\rm Logarithmic:} y = b_{\rm LOG} + a_{\rm LOG}\log(e_{\rm n})$'
        else:
            math_str = r'${\rm Logarithmic:} e_{\rm n} = \exp\left( \frac{y - b_{\rm LOG}}{a_{\rm LOG}}\right) $'

    elif FDCPar['model'] == 'power':
        min_vals = np.array([0, 0])
        max_vals = np.array([np.inf, np.inf])
        units = ['(-)', '(L/T)']
        if FDCPar['form'] == 'y':
            math_str = r'${\rm Power:} y = b_{\rm PW}e_{\rm n}^{- a_{\rm PW} }$'
        else:
            math_str = r'${\rm Power:} e_{\rm n} = \left(\frac{1}{b_{\rm PW}}y\right)^{(-1/a_{\rm PW})} $'

    elif FDCPar['model'] == 'quimpo':
        min_vals = np.array([0, -np.inf])
        max_vals = np.array([np.inf, np.inf])
        units = ['(L/T)', '(-)']
        if FDCPar['form'] == 'y':
            math_str = r'${\rm Quimpo:} y = a_{\rm Q} \exp\left(-b_{\rm Q}e_{\rm n}\right)$'
        else:
            math_str = r'${\rm Quimpo:} e_{\rm n} = - \frac{1}{b_{\rm Q}}\log\left(\frac{1}{a_{\rm Q}}y\right) $'

    elif FDCPar['model'] == 'viola':
        min_vals = np.array([0, 0])
        max_vals = np.array([np.inf, np.inf])
        units = ['(L/T)', '(-)']
        if FDCPar['form'] == 'y':
            math_str = r'${\rm Viola:} y = a_{\rm V} \left(\frac{1-e_{\rm n}}{e_{\rm n}}\right)^{b_{\rm V}}$'
        else:
            math_str = r'${\rm Viola:} e_{\rm n} = \left[\left(\frac{y}{a_{\rm V}}\right)^{(1/b_{\rm V})} + 1\right]^{-1} $'

    elif FDCPar['model'] == 'genuchten-2':
        min_vals = np.array([0, 0])
        max_vals = np.array([np.inf, np.inf])
        units = ['(T/L)', '(-)']
        if FDCPar['form'] == 'y':
            math_str = r'${\rm Genuchten-2:} y = \frac{1}{a_{\rm VG}} \left[e_{\rm n}^{\bigl(-b_{\rm VG}/(1-b_{\rm VG})\bigr)} - 1\right]^{(1/b_{\rm VG})}$'
        else:
            math_str = r'${\rm Genuchten-2:} e_{\rm n} = \left[1 + (a_{\rm VG}{y})^{b_{\rm VG}}\right]^{(1/b_{\rm VG}-1)}$'

    elif FDCPar['model'] == 'kosugi-2':
        min_vals = np.array([0, 0])
        max_vals = np.array([np.inf, np.inf])
        units = ['(L/T)', '(-)']
        if FDCPar['form'] == 'y':
            math_str = r'${\rm Kosugi-2:} y = a_{\rm K}\exp\left(\sqrt{2}b_{\rm K}{\rm erfc}^{-1}(2e_{\rm n})\right)$'
        else:
            math_str = r'${\rm Kosugi-2:} e_{\rm n} = \frac{1}{2}{\rm erfc}\left[\frac{1}{\sqrt{2} b_{\rm K}}\log\left(\frac{y}{a_{\rm K}}\right)\right]$'

    elif FDCPar['model'] == 'lognormal-3':
        min_vals = np.array([-np.inf, 0, -np.inf])
        max_vals = np.array([np.inf, np.inf, np.inf])
        units = ['(L/T)', '(L/T)', '(L/T)']
        if FDCPar['form'] == 'y':
            math_str = r'${\rm Lognormal-3:} y = c_{\rm LN} + \exp\left[ a_{\rm LN} - \sqrt{2}b_{\rm LN}{\rm erfc}^{-1}\bigl( 2(1-e_{\rm n})\bigr)\right]$'
        else:
#           math_str = r'${\rm Lognormal-3:} e_{\rm n} = \left\{\begin{array}{ll} 1 - \frac{1}{2}{\rm erfc}\left[\frac{a_{\rm LN} - \log(y - c_{\rm LN})}{\sqrt{2} b_{\rm LN}}\right] \hspace{0.5cm} {\rm if} \quad y > c_{\rm LN} \\ 1\hphantom{- \frac{1}{2}{\rm erfc}\left[\frac{a_{\rm LN} - \log(y - c_{\rm LN})}{\sqrt{2} b_{\rm LN}}\right]} \hspace{0.7cm} {\rm if} \quad y \leq c_{\rm LN} \end{array} \right.$'
           math_str = r'${\rm Lognormal-3:} e_{\rm n} = 1 - \frac{1}{2}{\rm erfc}\left[\frac{a_{\rm LN} - \log(y - c_{\rm LN})}{\sqrt{2} b_{\rm LN}}\right] \text{ if } y > c_{\rm LN} \text{ and } 1 \text{ if } y \leq c_{\rm LN}$'

    elif FDCPar['model'] == 'pareto':
        min_vals = np.array([-np.inf, 0, -np.inf])
        max_vals = np.array([np.inf, np.inf, np.inf])
        units = ['(L/T)', '(L/T)', '(-)']
        if FDCPar['form'] == 'y':
            math_str = r'${\rm Pareto:} y = a_{\rm GP} + \frac{b_{\rm GP}}{c_{\rm GP} } \left[e_{\rm n}^{-c_{\rm GP}} - 1 \right]$'
        else:
#            math_str = r'${\rm Pareto:} e_{\rm n} = 1 - \left\{\begin{array}{ll} \left[1+ \frac{c_{\rm GP}(y - a_{\rm GP})}{b_{\rm GP}}\right]^{(-1/c_{\rm GP})} & {\rm if} \quad c_{\rm GP} \neq 0 \\ \exp\left(\frac{a_{\rm GP} - y}{b_{\rm GP}}\right) & {\rm if} \quad c_{\rm GP} = 0 \end{array} \right.$'
            math_str = r'${\rm Pareto:} e_{\rm n} = 1 - \left[1+ \frac{c_{\rm GP}(y - a_{\rm GP})}{b_{\rm GP}}\right]^{(-1/c_{\rm GP})} \text{ if } c_{\rm GP} \neq 0 \text{ and } \exp\left(\frac{a_{\rm GP} - y}{b_{\rm GP}}\right) \text{ if } c_{\rm GP} = 0$'

    elif FDCPar['model'] == 'gev':
        min_vals = np.array([-np.inf, 0, -np.inf])
        max_vals = np.array([np.inf, np.inf, np.inf])
        units = ['(L/T)', '(L/T)', '(-)']
        if FDCPar['form'] == 'y':
            math_str = r'${\rm GEV:} y = a_{\rm GEV} + \frac{b_{\rm GEV}}{c_{\rm GEV}}\left[\log\left(\frac{1}{(1 - e_{\rm n})}\right)^{-c_{\rm GEV}} - 1\right]$'
        else:
#            math_str = r'${\rm GEV:} e_{\rm n} = \left\{\begin{array}{ll} 1 - \exp\left\{-\left[1 + c_{\rm GEV} \left(\frac{y - a_{\rm GEV}}{b_{\rm GEV}}\right)\right]^{(-1/c_{\rm GEV})}\right\} & {\rm if} \quad c_{\rm GEV} \neq 0 \\ 1 - \exp \left[ -\exp \left(\frac{a_{\rm GEV} - y}{b_{\rm GEV}}\right)\right] & {\rm if} \quad c_{\rm GEV} = 0 \end{array} \right.$'
            math_str = r'${\rm GEV:} e_{\rm n} = 1 - \exp\left\{-\left[1 + c_{\rm GEV} \left(\frac{y - a_{\rm GEV}}{b_{\rm GEV}}\right)\right]^{(-1/c_{\rm GEV})}\right\} \text{ if } c_{\rm GEV} \neq 0 \text{ and } 1 - \exp \left[ -\exp \left(\frac{a_{\rm GEV} - y}{b_{\rm GEV}}\right)\right] \text{ if } c_{\rm GEV} = 0$'

    elif FDCPar['model'] == 'franchini':
        min_vals = np.array([0, -np.inf, 0])
        max_vals = np.array([np.inf, np.inf, np.inf])
        units = ['(L/T)', '(L/T)', '(-)']
        if FDCPar['form'] == 'y':
            math_str = r'${\rm Franchini:} y = b_{\rm FS} + a_{\rm FS}(1 - e_{\rm n})^{c_{\rm FS}}$'
        else:
            math_str = r'${\rm Franchini:} e_{\rm n} = 1 - \left(\frac{y - b_{\rm FS}}{a_{\rm FS}}\right)^{(1/c_{\rm FS})}$'

    elif FDCPar['model'] == 'genuchten-3':
        min_vals = np.array([0, 1, 0])
        max_vals = np.array([np.inf, np.inf, 1])
        units = ['(T/L)', '(-)', '(-)']
        if FDCPar['form'] == 'y':
            math_str = r'${\rm Genuchten-3:} y = \frac{1}{a_{\rm VG}} \left[e_{\rm n}^{(-1/c_{\rm VG})} - 1 \right]^{(1/b_{\rm VG})}$'
        else:
            math_str = r'${\rm Genuchten-3:} e_{\rm n} = \left[ 1 + (a_{\rm VG} y)^{b_{\rm VG}} \right]^{-c_{\rm VG}}$'

    elif FDCPar['model'] == 'kosugi-3':
        min_vals = np.array([0, 0, -np.inf])
        max_vals = np.array([np.inf, np.inf, np.inf])
        units = ['(L/T)', '(-)', '(L/T)']
        if FDCPar['form'] == 'y':
            math_str = r'${\rm Kosugi-3:} y = c_{\rm K} + (a_{\rm K} - c_{\rm K})\exp\left[\sqrt{2} b_{\rm K} \text{erfc}^{-1}(2e_{\rm n})\right]$'
        else:
            math_str = r'${\rm Kosugi-3:} e_{\rm n} = \frac{1}{2} \, \text{erfc}\left[\frac{1}{\sqrt{2} b_{\rm K}}\log\left(\frac{y - c_{\rm K}}{a_{\rm K} - c_{\rm K}}\right)\right] \text{ if } y > c_{\rm K} \text{ and } e_{\rm n} = 1 \text{ otherwise}$'

    elif FDCPar['model'] == 'weibull':
        min_vals = np.array([0, 0])
        max_vals = np.array([np.inf, np.inf])
        units = ['(L/T)', '(L/T)']
        math_str = ''  # No mathematical expression provided in the original MATLAB code

    elif FDCPar['model'] == 'non_central_F':
        min_vals = np.array([0, 0, -np.inf])
        max_vals = np.array([np.inf, np.inf, np.inf])
        units = ['(-)', '(L/T)', '(L/T)']
        math_str = ''  # No mathematical expression provided in the original MATLAB code

    elif FDCPar['model'] == 'exponential':
        min_vals = np.array([0, 0])
        max_vals = np.array([np.inf, np.inf])
        units = ['(-)', '(L/T)']
        math_str = ''  # No mathematical expression provided in the original MATLAB code

    elif FDCPar['model'] == 'gamma':
        min_vals = np.array([0, 0])
        max_vals = np.array([np.inf, np.inf])
        units = ['(-)', '(L/T)']   # b enters as scale parameter
        math_str = ''  # No mathematical expression provided in the original MATLAB code

    elif FDCPar['model'] == 'rayleigh':
        min_vals = np.array([0])
        max_vals = np.array([np.inf])
        units = ['(L/T)']
        math_str = ''  # No mathematical expression provided in the original MATLAB code

    else:
        raise ValueError("Unknown model")

    # Store parameter range information
    Par_info = {'min': min_vals.copy(), 'max': max_vals.copy()}
    # Copy as entries of Par_info will change if min_vals/max_vals change
    
    # Ensure default settings for options if not specified
    if 'TolX' not in options:
        options['TolX'] = 1e-2
    if 'TolFun' not in options:
        options['TolFun'] = 1e-3
    if 'MaxFunEvals' not in options:
        options['MaxFunEvals'] = 1e4
    
    # Algorithm-specific settings for LM/SP
    if method.lower() in ['sp', 'lm']:
        if 'N' not in options:
            options['N'] = int(5)
    
    # Algorithm-specific settings for DE/CMAES
    if method.lower() in ['de', 'cma']:
        if 'P' not in options:
            options['P'] = int(25)
    
    # Verify DE-specific settings
    if method.lower() == 'de':
        if 'CR' not in options:
            options['CR'] = 0.8
    
    # Set LSQ options for least squares optimization methods
    if method.lower() in ['de', 'cmaes']:
        lsq_options = None
    elif method.lower() == 'lm':
        lsq_options = {
            'disp': 0,
            'xtol': options['TolX'],
            'ftol': options['TolFun'],
            'max_nfev': int(options['MaxFunEvals'])
            }
    else: # for Nelder-Mead simplex algorithm
        lsq_options = {
            'disp': 'False',
            'xatol': options['TolX'],
            'fatol': options['TolFun'],
            'maxfev': int(options['MaxFunEvals'])
            }

    # Handle specific method options
    if method.lower() not in ['de', 'cma']:
        if 'P' in options:
            del options['P']
    if method.lower() not in ['lm', 'sp']:
        if 'N' in options:
            del options['N']
    
    # Remove DE-specific fields if method is not 'de'
    if method.lower() != 'de':
        for field in ['CR', 'F']:
            if field in options:
                del options[field]
    
    if 'print' not in options:
        options['print'] = 'yes'
    
    if 'type' not in options:
        options['type'] = 'daily'
    
    # Next we set the ranges for the initial parameter values
    for i in range(len(min_vals)):
        if np.isinf(min_vals[i]):
            min_vals[i] = -1    # 
        if np.isinf(max_vals[i]):
            max_vals[i] = 2     # used to be 1 but leads to problems VG-3 as b>1 for c>0
   
    # -> still have from -inf to inf if allowed
    # -> can pose problems with boundary enforcement --> reflection so that it only uses 1 bound

    # Generate samples based on method
    if method.lower() == 'de':  # Draw initial population for DE
        X = LH_sampling(min_vals, max_vals, options['P'])
    elif method.lower() == 'lm':  # Draw initial population for LM
        X = LH_sampling(min_vals, max_vals, options['N'])
    elif method.lower() == 'sp':  # Draw multistart points for SP
        X = LH_sampling(min_vals, max_vals, options['N'])
    else:  # CMAES creates its own population
        X = []

    # Print summaries to screen
    print('---------------- Summary of fields FDCPar -----------------')
    for key, value in FDCPar.items():
        print(f"{key:<7}: {value:<6}")
    print('-----------------------------------------------------------')

    print('--------------- Summary of fields options -----------------')
    for key, value in options.items():
        print(f"{key:<11}: {value:<6}")
    print('-----------------------------------------------------------')
     
    return X, FDCPar, Par_info, options, lsq_options, n_crash, units, math_str


def FDCFIT_functions(x, FDCPar, e, y, returnarg):
    # ####################################################################### #
    #                                                                         #
    # This function evaluates the exceedance probability, e_s, or streamflow, #
    # y_s, using a suite of different parametric expressions of the flow      #
    # duration curve. Return argument is set by the user by the last input    #
    # argument, 1: L2 objective function, 2: residual vector, e-e_s or y-y_s, #
    # and 3: L1 objective function                                            # 
    #                                                                         #
    # SYNOPSIS: varargout = FDCFIT_functions(x,FDCPar,e,y,returnarg)          #
    #  where                                                                  #
    #   x         [input] 1xd vector parameter values parametric expression   #
    #   FDCPar    [input] Structure with settings for FDCFIT                  #
    #   e         [input] nx1 vector of empirical exceedance probabilities    #
    #   y         [input] nx1 vector of empirical (= measured) discharges     #
    #   returnarg [input] Integer of desired return argument                  #
    #    = 1            L2 objective function                                 #
    #    = 2            Residual vector of exceedance or discharge values     #
    #    = 3            L1 objective function                                 #
    #   varargout [outpt] L2/L1 objective function or residual vector         #
    #                                                                         #
    # (c) Written by Jasper A. Vrugt, Dec. 2014                               #
    # University of California Irvine                                         #
    #                                                                         #
    # ####################################################################### #

    # Extract parameters
    a = x[0]
    b = x[1]

    # Switch between different models based on FDCPar['model']
    # Two-parameter formulations
    if FDCPar['model'] == 'lognormal-2':
        # Lognormal distribution with 2 parameters
        if FDCPar['form'] == 'e':  # If exceedance probability space (e)
            e_s = 1 - 0.5 * erfc(-(np.log(y) - a) / (np.sqrt(2) * b))
        else:  # If discharge space (y)
            y_s = np.exp(a - np.sqrt(2) * b * erfcinv(2 * (1 - e)))
    
    elif FDCPar['model'] == 'gumbel':
        # Gumbel distribution with 2 parameters
        if FDCPar['form'] == 'e':
            e_s = 1 - np.exp(-np.exp(-(y - a) / b))
        else:
            y_s = a - b * np.log(-np.log(1 - e))
    
    elif FDCPar['model'] == 'logistic':
        # Logistic distribution with 2 parameters
        if FDCPar['form'] == 'e':
            e_s = 1 - 1 / (1 + np.exp(-(y - a) / b))
        else:
            y_s = a - b * np.log(1 / (1 - e) - 1)
    
    elif FDCPar['model'] == 'logarithmic':
        # Logarithmic distribution with 2 parameters
        if FDCPar['form'] == 'e':
            e_s = np.exp((y - b) / a)
        else:
            y_s = b + a * np.log(e)
    
    elif FDCPar['model'] == 'power':
        # Power distribution with 2 parameters
        if FDCPar['form'] == 'e':
            e_s = (y / b) ** (-1 / a)
        else:
            y_s = b * (e) ** (-a)
    
    elif FDCPar['model'] == 'quimpo':
        # Quimpo distribution with 2 parameters
        if FDCPar['form'] == 'e':
            e_s = -1 / b * np.log(y / a)
        else:
            y_s = a * np.exp(-b * e)
    
    elif FDCPar['model'] == 'viola':
        # Viola distribution with 2 parameters
        if FDCPar['form'] == 'e':
            e_s = 1 / ((y / a) ** (1 / b) + 1)
        else:
            y_s = a * (1 / e - 1) ** b
    
    elif FDCPar['model'] == 'genuchten-2':
        # Van Genuchten distribution with 2 parameters
        c = 1 - 1 / b
        # Protection for out of bound - only this and VG-3 expression
        if any(param < 0 for param in [a, b, c]):
            e_s, y_s = np.zeros_like(e), np.zeros_like(y)  # Initialize arrays with zeros
        else:
            if FDCPar['form'] == 'e':
                e_s = (1 + (a * y) ** b) ** (-c)
            else:
                y_s = (1 / a) * (e ** (-1 / c) - 1) ** (1 / b)
    
    elif FDCPar['model'] == 'kosugi-2':
        # Kosugi distribution with 2 parameters
        c = 0
        if FDCPar['form'] == 'e':
            e_s = np.ones_like(y)
            dummy = (y - c) / (a - c)
            id_ = dummy > 0
            e_s[id_] = 1 / 2 * erfc(np.log((y[id_] - c) / (a - c)) / (np.sqrt(2) * b))
        else:
            y_s = c + (a - c) * np.exp(erfcinv(2 * e) * np.sqrt(2) * b)

    # Three-parameter formulations
    elif FDCPar['model'] == 'lognormal-3':
        # Lognormal distribution with 3 parameters
        c = x[2]
        if FDCPar['form'] == 'e':
            e_s = np.ones_like(y)
            id = np.where(y > c)
            e_s[id] = 1 - 0.5 * erfc(-(np.log(y[id] - c) - a) / (np.sqrt(2) * b))
        else:
            y_s = c + np.exp(a - np.sqrt(2) * b * erfcinv(2 * (1 - e)))
        
    elif FDCPar['model'] == 'pareto':
        # Generalized Pareto distribution with 3 parameters
        c = x[2]
        if FDCPar['form'] == 'e':
            theta = np.maximum((y - a) / b, 0)
            cdf_p = 1 - (1 + c * theta) ** (-1 / c)
            if c < 0:
                cdf_p[theta > -1 / c] = 1
            e_s = 1 - cdf_p
        else:
            y_s = a + b / c * (e ** (-c) - 1)
            y_s[e > 1] = 0
    
    elif FDCPar['model'] == 'gev':
        # Generalized Extreme Value distribution with 3 parameters
        c = x[2]
        if FDCPar['form'] == 'e':
            if c > 0:
                y_check = a - b / c
                id = np.where(y < y_check)
                value = 0
            else:
                y_check = a + b / (-c)
                id = np.where(y > y_check)
                value = 1
            cdf_gev = np.exp(- (1 + c * ((y - a) / b)) ** (-1 / c))
            cdf_gev[id] = value
            e_s = 1 - cdf_gev
        else:
            if c > 0:
                y_check = a - b / c
                y_s = a + b / c * ((-np.log(1 - e)) ** (-c) - 1)
                y_s[y < y_check] = a - b / c
            else:
                y_check = a + b / (-c)
                y_s = a + b / c * ((-np.log(1 - e)) ** (-c) - 1)
                y_s[y > y_check] = np.inf

    elif FDCPar['model'] == 'franchini':
        # Franchini expression with 3 parameters
        c = x[2]
        if FDCPar['form'] == 'e':
            e_s = np.ones_like(y)
            dummy = (y - b) / a
            id = np.where(dummy > 0)
            e_s[id] = 1 - dummy[id] ** (1 / c)
        else:
            y_s = b + a * (1 - e) ** c

    elif FDCPar['model'] == 'genuchten-3':
        # van Genuchten expression with 3 parameters
        c = x[2]
        if any(param < 0 for param in [a, b, c]):
            e_s, y_s = np.zeros_like(e), np.zeros_like(y)  # Initialize arrays with zeros
        else:
            if FDCPar['form'] == 'e':
                e_s = (1 + (a * y) ** b) ** (-c)
            else:
                y_s = (1 / a) * (e ** (-1 / c) - 1) ** (1 / b)

    elif FDCPar['model'] == 'kosugi-3':
        # Kosugi expression with 3 parameters
        c = x[2]
        if FDCPar['form'] == 'e':
            e_s = np.ones_like(y)
            dummy = (y - c) / (a - c)
            id = np.where(dummy > 0)
            e_s[id] = 0.5 * erfc(np.log(dummy[id]) / (np.sqrt(2) * b))
        else:
            y_s = c + (a - c) * np.exp(erfcinv(2 * e) * np.sqrt(2) * b)

    # Alternative models not used in original paper
    # Three parameter formulations    
    elif FDCPar['model'] == 'weibull':
        # Weibull distribution with 2 parameters
        if FDCPar['form'] == 'e':
            e_s = np.exp(-(y / a) ** b)
        else:
            y_s = -a * (np.log(e)) ** (1 / b)

    elif FDCPar['model'] == 'non_central_F':
        # Non-central-F distribution with 3 parameters
        c = x[2]
        if FDCPar['form'] == 'e':
            e_s = 1 - stats.ncf.cdf(y,a,b,c) # ncfd.cdf(y, a, b, c)
        else:
            raise ValueError("Not available in 'y' form")

    # Alternative models not used in original paper
    # Two-parameter formulations
    elif FDCPar['model'] == 'exponential':
        # Exponential distribution with 2 parameters
        if FDCPar['form'] == 'e':
            e_s = (1 / a) * np.log(y / b)
        else:
            y_s = b * np.exp(a * e)
    
    elif FDCPar['model'] == 'gamma':
        # Gamma distribution with 2 parameters
        if FDCPar['form'] == 'e':
            e_s = 1 - gamma.cdf(y, a, scale=b)
        else:
            y_s = gamma.ppf(1 - e, a, scale=b)

    elif FDCPar['model'] == 'rayleigh':
        # Rayleigh distribution with 1 parameter!
        if FDCPar['form'] == 'e':
            e_s = np.exp(- 0.5 * (y ** 2 / a ** 2))
        else:
            y_s = np.sqrt(-2 * a ** 2 * np.log(e))

    else:
        raise ValueError("Unknown model")

    # Residual computation
    if FDCPar['form'] == 'e':
        res = e - e_s
    else:
        res = y - y_s
    # Return based on returnarg value
    if returnarg == 0:
        if FDCPar['form'] == 'e':
            return e_s
        else:
            return y_s    
    elif returnarg == 1:
        return np.dot(res, res)  # L2 objective function (sum of squares)
    elif returnarg == 2:
        return res  # Residual vector
    elif returnarg == 3:
        return np.linalg.norm(res, 1)  # L1 objective function (sum of absolute values)
    else:
        raise ValueError("FDCFIT_functions: Unknown return argument")


def FDCFIT_end(FDCPar, MAP, RMSE_MAP, it_MAP):
    # ####################################################################### #
    #                                                                         #
    # Finalize return arguments, file writing, and setup calculation          #
    #                                                                         #
    # SYNOPSIS: [map,RMSE_map,it_map,str] = FDCFIT_end(FDCPar,MAP, ...        #
    #               RMSE_MAP,it_MAP)                                          #
    #  where                                                                  #
    #   FDCPar    [input] Structure with settings for FDCFIT                  #
    #   MAP       [input] Nxd matrix of N optimized parameter vectors         #
    #   RMSE_MAP  [input] Nx1 vector of RMSE of N parameter vectors           #
    #   it_MAP    [input] # iterations each of N multi-start trials           # 
    #   map       [outpt] 1xd vector with overall best (map) parameter values #
    #   RMSE_map  [outpt] scalar with RMSE of map parameter values            #
    #   it_map    [outpt] # iterations of map                                 #
    #   str       [outpt] string with parameter names                         #
    #                                                                         #
    # (c) Written by Jasper A. Vrugt, Dec. 2014                               #
    # University of California Irvine                                         #
    #                                                                         #
    # ####################################################################### #
 
    # Find the index of the minimum RMSE value
    id_best = np.nanargmin(RMSE_MAP)

    # Extract the best parameter values and corresponding RMSE
    map = MAP[id_best, :FDCPar['d']]  # Extract the best parameter vector
    RMSE_map = RMSE_MAP[id_best]  # RMSE of the best parameters
    it_map = it_MAP[id_best]  # Number of iterations for the best parameters

    # Define the parameter names as strings
    str = ['a_', 'b_', 'c_']
    
    # Open an output file for warnings and append
    with open('warning_file.txt', 'a+') as fid:
        # Write final line to the warning file
        fid.write('----------- End of MODELAVG warning file ----------\n')
    
    # On Windows or macOS, prompt to edit the warning file
    if any(platform in sys.platform for platform in ["win", "darwin"]):
        # Open the warning file with the default text editor
        os.startfile('warning_file.txt')

    # Return the best parameters, RMSE, iterations, and parameter names
    return map, RMSE_map, it_map, str


def FDCFIT_plot(FDCPar, e, y, method, map, RMSE_map, str, units, math_str):
    # ####################################################################### #
    #                                                                         #
    # Plots the results of the flow duration curve formulation                #
    #                                                                         #
    # SYNOPSIS: FDCFIT_plot(FDCPar,e,y,method,map,RMSE_map,str,units, ...     #
    #               math_str)                                                 #
    #  where                                                                  #
    #   FDCPar    [input] Structure with settings for FDCFIT                  #
    #   e         [input] nx1 vector of empirical exceedance probabilities    #
    #   y         [input] nx1 vector of empirical (= measured) discharges     #
    #   method    [input] Name (string) of parameter estimation method        #
    #   map       [input] 1xd vector optmzed parameters parametric FDC expr.  #
    #   RMSE_map  [input] scalar with minimized objective function            #
    #   str       [input] string with model equation FDC expression in latex  #
    #   units     [input] string with model equation FDC expression in latex  #
    #   math_str  [input] string with model equation FDC expression in latex  #
    #                                                                         #
    # (c) Written by Jasper A. Vrugt, Dec. 2014                               #
    # University of California Irvine                                         #
    #                                                                         #
    # ####################################################################### #
    
     # Define name of program
    n_program = 'FDCFIT'
    
    # Define name of figures file
    file_name = f'{n_program}_figures.pdf'

    # Delete the warning file and open output file for results
    try:
        with open('FDCFIT_output.txt', 'w', encoding='utf-8') as fid:
            fid.write('------------------- FDCFIT output file -------------------\n')
            fid.write('\n')

            if FDCPar['form'] == 'e':
                fid.write(f"→ Results for model {FDCPar['model'].upper()} (Check Table 2 of manual)\n")
            else:
                fid.write(f"→ Results for model {FDCPar['model'].upper()} (Check Table 1 of manual)\n")

            if FDCPar['form'] == 'y':
                fid.write("→ Flow duration curve minimizes discharge residuals\n")
            else:
                fid.write("→ Flow duration curve minimizes exceedance probability residuals\n")
            
            fid.write(f"→ Parameter estimation method: {method.upper()} algorithm\n")
            fid.write('\n')
            fid.write("==================================\n")
            
            # Print parameter values to file
            for i in range(FDCPar['d']):
                param_name = f"{str[i]}{FDCPar['acronym']}"
                value = map[i]
                unit = units[i]
                fid.write(f" {param_name:<5}   {value:>10.5f}  {unit:<5}\n")

            fid.write('----------------------------------\n')
            if FDCPar['form'] == 'e':
                unit = '(-)'
                fid.write(f" RMSE    {RMSE_map:>10.5f}  {unit}   --> in exceedance probability space\n")
            else:
                unit = '(L/T)'
                fid.write(f" RMSE    {RMSE_map:>10.5f}  {unit}   --> in discharge space\n")    
            fid.write("==================================\n")
            fid.write("\n")
            fid.write('---------------- End of FDCFIT output file -------------------\n')
    
    except IOError as e:
        print(f"Error writing output file: {e}")
    
    # Open warning file if on Windows or Mac
    if sys.platform == "win32" or sys.platform == "darwin":
        os.startfile("FDCFIT_output.txt")

    # Determine screen size (using matplotlib to get screen dimensions)
    monitor = get_monitors()[0]
    screen_width = monitor.width
    screen_height = monitor.height
    x_mult = screen_width / 1920
    y_mult = screen_height / 1080
    t_mult = min(x_mult, y_mult)

    # Define fontsize for figures
    fontsize_xylabel = 18 * t_mult
    fontsize_axis = 16 * t_mult
    fontsize_legend = 14 * t_mult
    fontsize_text = 14 * t_mult
    fontsize_table = 16 * t_mult
    fontsize_titlepage = 20 * t_mult
    
    # ----------------------------------------------------------------------- %
    # Now plot empty figure for PDF file
    # ----------------------------------------------------------------------- %

    plt.rcParams['text.usetex'] = False
    with PdfPages(file_name) as pdf:
        # Create an empty figure for the PDF
        plt.figure(figsize=(12, 8))
        plt.plot([], [], 'ro')
        plt.axis([0, 1, 0, 1])
        plt.gca().set_facecolor('w')
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])
        plt.text(0.3 * x_mult, 0.6 * y_mult, r"${\rm Visual \; results \; of \; FDCFIT \; toolbox}$", fontsize=fontsize_titlepage)
        pdf.savefig()
        plt.close()
       
        # Determine the scaling for the plot
        Y_max = np.floor(1.02 * np.max(y))
        Y_min = max(1e-4, np.min(y))
        X_min, X_max = -0.04, 1.04
    
        # Now plot the measured FDC
        plt.figure(figsize=(15, 10))
 
        # Handle second plot with linear y-scale
        ax = plt.subplot(1,2,1)
        ax.plot(e, y, 'ro', markersize=8)  # Measured data in red
        if FDCPar['form'] == 'e':
            e_s = FDCFIT_functions(map, FDCPar, e, y, 0)
            ax.plot(e_s, y, 'b-', linewidth=2)  # Model fit in blue
        else:
            y_s = FDCFIT_functions(map, FDCPar, e, y, 0)
            ax.plot(e, y_s, 'b-', linewidth=2)  # Model fit in blue

        ax.set_xticks(np.arange(0, 1.1, 0.1))
        ax.set_xticklabels([f'{x:.1f}' for x in np.arange(0, 1.1, 0.1)])
        ax.set_xlim([X_min, X_max])
        ax.set_ylim([Y_min, Y_max])
        ax.tick_params(axis='both', labelsize=fontsize_axis)
        # Add RMSE text
        if FDCPar['form'] == 'e':
            ax.text(0.12, 0.8, r"$\text{RMSE}:$"f'{RMSE_map:.5f}'r"$\text{ (-)}$", transform=ax.transAxes, fontsize=fontsize_text) # fontsize_legend, color='black', verticalalignment='center')
        else:
            ax.text(0.12, 0.8, r"$\text{RMSE}:$"f'{RMSE_map:.5f}'r"$\text{ (L/T)}$", transform=ax.transAxes, fontsize=fontsize_text) # fontsize_legend, color='black', verticalalignment='center')
        ax.set_xlabel(r'${\rm Normalized \, exceedance \, probability}, \, e_\text{n}\;(-)$', fontsize=fontsize_xylabel)
        ax.set_ylabel(r'${\rm Streamflow}, \, y \,\, {\rm (L/T)}$', fontsize=fontsize_xylabel)                     
        ax.text(0., 1.02, r"${\rm (a) \, Linear \, scale}$", transform=ax.transAxes, fontsize=fontsize_titlepage) 
 
        # Handle second plot with logarithmic y-scale
        ax = plt.subplot(1, 2, 2)
        ax.semilogy(e, y, 'ro', markersize=8)  # Measured data in red

        if FDCPar['form'] == 'e':
            ax.semilogy(e_s, y, 'b-', linewidth=2)  # Model fit in blue
        else:
            ax.semilogy(e, y_s, 'b-', linewidth=2)  # Model fit in blue
        ax.set_xticks(np.arange(0, 1.1, 0.1))
        ax.set_xticklabels([f'{x:.1f}' for x in np.arange(0, 1.1, 0.1)])
        ax.set_xlim([X_min, X_max])
        ax.set_ylim([Y_min, Y_max])
        ax.tick_params(axis='both', labelsize=fontsize_axis)
        # Add RMSE to the plot
        if FDCPar['form'] == 'e':
            ax.text(0.32, 0.8, r"$\text{RMSE}:$"f'{RMSE_map:>9.5f}'r"$\text{ (-)}$", transform=ax.transAxes, fontsize=fontsize_text) # fontsize_legend, color='black', verticalalignment='center')
        else:
            ax.text(0.32, 0.8, r"$\text{RMSE}:$"f'{RMSE_map:>9.5f}'r"$\text{ (L/T)}$", transform=ax.transAxes, fontsize=fontsize_text) # fontsize_legend, color='black', verticalalignment='center')
        # Add parameter values to the plot
        for i in range(FDCPar['d']):
            param_name = f"{str[i]}{FDCPar['acronym']}"
            value = map[i]
            unit = units[i]
            param_str = f"{param_name:>5}: {value:>10.5f} {unit}"     
            # Display the text on the plot
            ax.text(0.32, 0.77 - i*0.03,param_str, transform=ax.transAxes, fontsize=fontsize_text) # fontsize_legend, color='black', verticalalignment='center')
        
        ax.set_xlabel(r'${\rm Normalized \, exceedance \, probability}, \, e_\text{n}\;(-)$', fontsize=fontsize_xylabel)
        ax.set_ylabel(r'${\rm Streamflow}, \, y \,\, {\rm (L/T)}$', fontsize=fontsize_xylabel)                     
        ax.text(0., 1.02, r"${\rm (b) \, Logarithmic \, scale}$", transform=ax.transAxes, fontsize=fontsize_titlepage) # fontsize_legend, color='black', verticalalignment='center')
        # Show legend
        ax.legend(['Measured data', math_str],fontsize=fontsize_legend, loc='upper right')
        pdf.savefig()
 
    plt.close()

    # Open the final PDFs
    os.startfile(file_name)

    # Print completion message
    print('FDCFIT: done')


def de_code(X, FDCPar, e, y, Par_info, options):
    # ####################################################################### #
    #                                                                         #
    # Differential Evolution Code                                             #
    #                                                                         #
    # SYNOPSIS: [X,RMSE,it] = de_code(X,FDCPar,e,y,Par_info,options)          #
    #                                                                         #
    # (c) Written by Jasper A. Vrugt, Dec. 2014                               #
    # University of California Irvine                                         #
    #                                                                         #
    # ####################################################################### #
    
    # Replicate parameter boundaries
    Par_info['min'] = np.tile(Par_info['min'], (options['P'], 1))
    Par_info['max'] = np.tile(Par_info['max'], (options['P'], 1))

    # Start counter
    it = options['P']
    tr = int(np.ceil(0.8 * options['P']))
    
    # Calculate initial objective function (OF)
    OF = np.full(options['P'], np.nan)
    for i in range(options['P']):
        OF[i] = FDCFIT_functions(X[i, :FDCPar['d']], FDCPar, e, y, 1)
    
    # Sort the population based on the objective function (OF)
    id_s = np.argsort(OF)
    F = 0.6  # Differential weight
    G = 0.4  # Crossover weight
    count = 0
    converged = False
    
    # Loop until converged or max evaluations reached
    while not converged and it < options['MaxFunEvals']:
        # Draw values to generate a, b, and c from
#        draw = np.sort(np.random.rand(options['P'] - 1, options['P']), axis=0)
        draw = np.argsort(np.random.rand(options['P'] - 1, options['P']), axis=0)
        # Create a, b, and c
        r = draw[:3, :].T  # Shape should be P x 3
        Z = X[:options['P'], :FDCPar['d']] + F * (X[r[:, 0], :FDCPar['d']] - X[r[:, 1], :FDCPar['d']]) + \
            G * (X[id_s[0], np.newaxis] - X[:options['P'], :FDCPar['d']])
        # Apply crossover
        crossover_mask = np.random.rand(options['P'], FDCPar['d']) > options['CR']
        Z[crossover_mask] = X[crossover_mask]
        # Make sure that everything is in bounds --> reflect if out of bounds
        Z[Z < Par_info['min']] = 2 * Par_info['min'][Z < Par_info['min']] - Z[Z < Par_info['min']]
        Z[Z > Par_info['max']] = 2 * Par_info['max'][Z > Par_info['max']] - Z[Z > Par_info['max']]
        
        # Calculate OF for each child
        OFc = np.full(options['P'], np.nan)
        for i in range(options['P']):
            OFc[i] = FDCFIT_functions(Z[i, :FDCPar['d']], FDCPar, e, y, 1)
        
        # Accept or reject new solutions
        accept_mask = OFc < OF
        X[accept_mask, :FDCPar['d']] = Z[accept_mask, :FDCPar['d']]
        OF[accept_mask] = OFc[accept_mask]
        
        # Check for convergence
        OF_s = np.sort(OF)
        id_s = np.argsort(OF)
        ii = id_s[:tr]
        
        it += options['P']
        
        # Print progress
        if it > 1:
            print('\r' + f"FDCFIT CALCULATING: {100 * (it / options['MaxFunEvals']):3.2f}% done of "
                  f"maximum of {options['MaxFunEvals']} trials (= options.MaxFunEvals) with Differential Evolution",
                  end='')
        
        # Check for convergence
        if (np.max(np.abs(np.max(X[ii, :], axis=0) - np.min(X[ii, :], axis=0))) < options['TolX']) and \
           ((OF_s[tr] - OF_s[0]) < options['TolFun']):
            converged = True
            print(f"\nFDCFIT CALCULATING: Differential Evolution has converged after {it} generations.")
        
    # Calculate RMSE
    RMSE = np.sqrt(OF / FDCPar['n'])
    # Return the number of iterations (repeated for all individuals)
    it = np.full(options['P'], it)

    return X, RMSE, it


def FDCFIT_jac(FDCPar, x, y):
    # ####################################################################### #
    #                                                                         #
    # Calculates the Jacobian of the VG or Kosugi FDC formulations using      #
    # analytical functions of the partial parameter sensitivities (zero cost) #
    #                                                                         #
    # SYNOPSIS: J = FDCFIT_jac(FDCPar,x,Y)                                    #
    #  where                                                                  #
    #   FDCPar    [input] Structure with settings for FDCFIT                  #
    #    .d           # unknown parameters of parametric FDC expression       #
    #    .form        Forward or inverse formulation of the FDC               #
    #      = 'e'      Forward expression: e_s = f(y)                          #
    #      = 'y'      Backward/inverse expression: y_s = f(e)                 #
    #   x         [input] 1xd vector parameter values parametric expression   #
    #   y         [input] nx1 vector of empirical (= measured) discharges     #
    #   J         [outpt] nxd Jacobian matrix of the FDC expression           #
    #                                                                         #
    # (c) Written by Jasper A. Vrugt, Dec. 2014                               #
    # University of California Irvine                                         #
    #                                                                         #
    # ####################################################################### #

    y = np.asarray(y)
    
    # For Van Genuchten (VG) model
    if FDCPar['model_class'] == 'vg':
        if FDCPar['model_formulation'] == '2':  # Two-parameter Van Genuchten
            a = x[0]
            b = x[1]
            c = 1 - 1 / b
            J = np.zeros((len(y), 2))
            J[:, 0] = (c * b * y * (a * y) ** (b - 1)) / ((a * y) ** b + 1) ** (c + 1)
            J[:, 1] = (c * np.log(a * y) * (a * y) ** b) / ((a * y) ** b + 1) ** (c + 1)

        elif FDCPar['model_formulation'] == '3':  # Three-parameter Van Genuchten
            a = x[0]
            b = x[1]
            c = x[2]
            J = np.zeros((len(y), 3))
            J[:, 0] = (c * b * y * (a * y) ** (b - 1)) / ((a * y) ** b + 1) ** (c + 1)
            J[:, 1] = (c * np.log(a * y) * (a * y) ** b) / ((a * y) ** b + 1) ** (c + 1)
            J[:, 2] = np.log((a * y) ** b + 1) / ((a * y) ** b + 1) ** c

        elif FDCPar['model_formulation'] == '5':  # Five-parameter mixture (Durner 1990)
            a = x[0]
            b = x[1]
            m_1 = 1 - 1 / b
            c = x[2]
            d = x[3]
            m_2 = 1 - 1 / d
            w_1 = x[4]
            w_2 = 1 - w_1
            J = np.zeros((len(y), 5))
            J[:, 0] = (m_1 * b * y * w_1 * (a * y) ** (b - 1)) / ((a * y) ** b + 1) ** (m_1 + 1)
            J[:, 1] = (m_1 * w_1 * np.log(a * y) * (a * y) ** b) / ((a * y) ** b + 1) ** (m_1 + 1)
            J[:, 2] = (m_2 * d * y * w_2 * (c * y) ** (d - 1)) / ((c * y) ** d + 1) ** (m_2 + 1)
            J[:, 3] = (m_2 * w_2 * np.log(c * y) * (c * y) ** d) / ((c * y) ** d + 1) ** (m_2 + 1)
            J[:, 4] = 1 / ((a * y) ** b + 1) ** m_1

        else:
            raise ValueError("Unknown formulation for Van Genuchten")

    # For Kosugi model
    elif FDCPar['model_class'] == 'k':
        if FDCPar['model_formulation'] == '2':  # Two-parameter Kosugi
            a = x[0]
            b = x[1]
            J = np.zeros((len(y), 2))
            J[:, 0] = -2 ** (1 / 2) / (2 * np.pi ** (1 / 2) * a * b * np.exp(np.log(y / a) ** 2 / (2 * b ** 2)))
            J[:, 1] = -(2 ** (1 / 2) * np.log(y / a)) / (2 * np.pi ** (1 / 2) * b ** 2 * np.exp(np.log(y / a) ** 2 / (2 * b ** 2)))

        elif FDCPar['model_formulation'] == '3':  # Three-parameter Kosugi
            a = x[0]
            b = x[1]
            c = x[2]
            J = np.zeros((len(y), 3))
            dummy = (y - c) / (a - c)
            id_valid = dummy > 0

            J[id_valid, 0] = -2 ** (1 / 2) / (2 * np.pi ** (1 / 2) * b * np.exp(np.log(-(c - y[id_valid]) / (a - c)) ** 2 / (2 * b ** 2)) * (a - c))
            J[id_valid, 1] = -(2 ** (1 / 2) * np.log(-(c - y[id_valid]) / (a - c))) / (2 * np.pi ** (1 / 2) * b ** 2 * np.exp(np.log(-(c - y[id_valid]) / (a - c)) ** 2 / (2 * b ** 2)))
            J[id_valid, 2] = (2 ** (1 / 2) * ((c - y[id_valid]) / (a - c) ** 2 + 1 / (a - c)) * (a - c)) / (2 * np.pi ** (1 / 2) * b * np.exp(np.log(-(c - y[id_valid]) / (a - c)) ** 2 / (2 * b ** 2)) * (c - y[id_valid]))

        elif FDCPar['model_formulation'] == '5':  # Five-parameter mixture (Durner 1990)
            a = x[0]
            b = x[1]
            psi_c1 = 0  # Psi_c1 is 0 by default
            c = x[2]
            d = x[3]
            w_1 = x[4]
            w_2 = 1 - w_1
            J = np.zeros((len(y), 5))
            J[:, 0] = -(2 ** (1 / 2) * w_1) / (2 * np.pi ** (1 / 2) * b * np.exp(np.log((psi_c1 - y) / (a - psi_c1)) ** 2 / (2 * b ** 2)) * (a - psi_c1))
            J[:, 1] = -(2 ** (1 / 2) * w_1 * np.log((psi_c1 - y) / (a - psi_c1))) / (2 * np.pi ** (1 / 2) * b ** 2 * np.exp(np.log((psi_c1 - y) / (a - psi_c1)) ** 2 / (2 * b ** 2)))
            J[:, 2] = -(2 ** (1 / 2) * w_2) / (2 * np.pi ** (1 / 2) * d * np.exp(np.log((psi_c1 - y) / (c - psi_c1)) ** 2 / (2 * d ** 2)) * (c - psi_c1))
            J[:, 3] = -(2 ** (1 / 2) * w_2 * np.log((psi_c1 - y) / (c - psi_c1))) / (2 * np.pi ** (1 / 2) * d ** 2 * np.exp(np.log((psi_c1 - y) / (c - psi_c1)) ** 2 / (2 * d ** 2)))
            J[:, 4] = -erfc((2 ** (1 / 2) * np.log((psi_c1 - y) / (a - psi_c1))) / (2 * b)) / 2

        else:
            raise ValueError("Unknown formulation for Kosugi")

    else:
        raise ValueError("Unknown water retention function")

    return J


def cmaes_code(FDCPar, e, y, options):
    # ####################################################################### #
    #                                                                         #
    # (µ/µ_w,lambda)-CMA-ES algorithm                                         #
    #                                                                         #
    # SYNOPSIS: [X,RMSE,it] = cmaes_code(FDCPar,e,y,options)                  #
    #                                                                         #
    # (c) Written by Jasper A. Vrugt, Dec. 2014                               #
    # University of California Irvine                                         #
    #                                                                         #
    # ####################################################################### #
    
    # Define constants
    d = FDCPar['d']         # number of objective variables/problem dimension
    x_mean = np.zeros(d)    # initial mean of objective variables
    sigma = 0.3             # coordinate-wise standard deviation (step size)
    stopfitness = 1e-10     # stop if fitness < stopfitness (minimization)
    stopeval = 1e4 * d**2   # stop after stopeval number of function evaluations

    # Strategy parameter setting: Selection
    lambda_ = options['P']
    mu = lambda_ / 2  # parents/points for recombination
    w = np.log(mu + 1 / 2) - np.log(np.arange(1, mu))  # weighted recombination
    mu = int(math.floor(mu))
    w /= np.sum(w)  # normalize recombination weights array
    mueff = np.sum(w)**2 / np.sum(w**2)  # variance-effectiveness of sum w_i x_i

    # Strategy parameter setting: Adaptation
    cc = (4 + mueff / d) / (d + 4 + 2 * mueff / d)  # time constant for cumulation for C
    cs = (mueff + 2) / (d + mueff + 5)  # t-const for cumulation for sigma control
    c1 = 2 / ((d + 1.3)**2 + mueff)  # learning rate for rank-one update of C
    cmu = min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((d + 2)**2 + mueff))  # rank-mu update
    damps = 1 + 2 * max(0, np.sqrt((mueff - 1) / (d + 1)) - 1) + cs  # damping for sigma

    # Initialize dynamic strategy parameters and constants
    pc = np.zeros(d)  # evolution path for C
    ps = np.zeros(d)  # evolution path for sigma
    B = np.eye(d)  # B defines the coordinate system
    D = np.ones(d)  # diagonal D defines the scaling
    C = B @ np.diag(D**2) @ B.T  # covariance matrix C
    invsqrtC = B @ np.diag(D**-1) @ B.T  # C^-1/2
    eigeneval = 0  # track update of B and D
    chiN = d**0.5 * (1 - 1 / (4 * d) + 1 / (21 * d**2))  # expectation of ||N(0,I)||

    child = np.nan * np.ones((d, lambda_))
    OFc = np.nan * np.ones(lambda_)
    it = 0
    count = 0

    while it < stopeval and it < options['MaxFunEvals']:
        # Generate and evaluate lambda offspring
        for k in range(lambda_):
            # m + sigma * Normal(0, C)
            child[:, k] = x_mean + sigma * B @ (D * np.random.randn(d))
            OFc[k] = FDCFIT_functions(child[:, k], FDCPar, e, y, 1)
            it += 1  # Update iteration counter
        # Sort by fitness and compute weighted mean into xmean
        # OFc_sorted, id_s = np.sort(OFc), np.argsort(OFc)
        OFc, id_s = np.sort(OFc), np.argsort(OFc)
        x_old = x_mean
        x_mean = child[:, id_s[:mu]] @ w  # recombination, new mean value
        # Cumulation: Update evolution paths
        ps = (1 - cs) * ps + np.sqrt(cs * (2 - cs) * mueff) * invsqrtC @ (x_mean - x_old) / sigma
        hsig = np.linalg.norm(ps) / np.sqrt(1 - (1 - cs)**(2 * it / lambda_)) / chiN < 1.4 + 2 / (d + 1)
        pc = (1 - cc) * pc + hsig * np.sqrt(cc * (2 - cc) * mueff) * (x_mean - x_old) / sigma

        tmp = (1 / sigma) * (child[:, id_s[:mu]] - np.tile(x_old, (mu, 1)).T)
        C = (1 - c1 - cmu) * C + c1 * (np.outer(pc, pc) + (1 - hsig) * cc * (2 - cc) * C) + cmu * tmp @ np.diag(w) @ tmp.T

        # Adapt step size sigma
        sigma = sigma * np.exp((cs / damps) * (np.linalg.norm(ps) / chiN - 1))

        # Decompose C into B * diag(D^2) * B'
        if (it - eigeneval) > (lambda_ / (c1 + cmu) / d / 10):
            eigeneval = it
            C = np.triu(C) + np.triu(C, 1).T  # enforce symmetry
            #B, D = np.linalg.eigh(C)   # eigen decomposition
            D, B = np.linalg.eigh(C)    # other way around in Python compared to matlab
            D = np.sqrt(D)
            #D = np.diag(D)              # make D a diagonal matrix as in MATLAB
            #D = np.sqrt(np.diag(D))     # D is vector standard deviations now
            invsqrtC = B @ np.diag(D**-1) @ B.T

        # Print progress
        if it > 1:
            print("\rFDCFIT CALCULATING: {:.2f}% done of maximum of {} trials (= options.MaxFunEvals) with CMA-ES algorithm".format(100 * (it / options['MaxFunEvals']), options['MaxFunEvals']), end="")

        # Check stopping conditions
        if OFc[0] <= stopfitness or np.max(D) > 1e7 * np.min(D) or (np.max(np.max(child, axis=1) - np.min(child, axis=1))) < options['TolX'] or (OFc[options['P'] - 1] - OFc[0]) < options['TolFun']:
            print(f"\nFDCFIT CALCULATING: CMAES has converged after {it} generations")
            break

    # Return parameter values
    X = child.T
    RMSE = np.sqrt(OFc / FDCPar['n'])
    it = np.full(options['P'], it)

    return X, RMSE, it


def LH_sampling(mn, mx, N):
    ## ################################################################################## ##
    ## Latin hypercube sampling of initial chain states DREAM Package                     ##
    ##                                                                                    ##
    ## SYNOPSIS: R = LH_sampling(mn,mx,N)                                                 ##
    ##  where                                                                             ##
    ##   mn        [input] REQUIRED: 1 x d vector of lower bound values                   ##
    ##   mx        [input] REQUIRED: 1 x d vector of upper bound values                   ##
    ##   N         [input] REQUIRED: # of Latin hypercube samples                         ##
    ##   r         [outpt] Nxd matrix of Latin hypercube samples                          ##
    ##                                                                                    ##
    ## Implementation based on the following paper                                        ##
    ##  Minasny, B., and A.B. McBratney (2006), A conditioned Latin hypercube method for  ##
    ##      sampling in the presence of ancillary information, Computers & Geosciences,   ##
    ##      Vol. 32 (9), pp. 1378-138                                                     ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Feb 2007                                             ##
    ## Los Alamos National Laboratory                                                     ##
    ##                                                                                    ##
    ## ################################################################################## ##

    d = len(mn)                                         # Number of parameters (dimensions)
    rng = np.array(mx) - np.array(mn)                   # 1 x d vector of parameter ranges
    y = np.random.rand(N, d)                            # N x d matrix with uniform random values
    # Python: important change below so that X in bound
    # as list is from 0 - N-1 rather than 1 to N
    # id_ = np.argsort(np.random.rand(N, d), axis=0)    # Sorting random values to avoid duplicates
    id_ = 1 + np.argsort(np.random.rand(N, d), axis=0)  # Random sort (1:N without replacement)
    M = (id_ - y) / N                                   # Multiplier matrix (introduces randomness)
    # Create the stratified LH samples
    R = np.add(mn, np.multiply(M, rng))                 # N x d matrix of stratified LH samples
    
    return R
