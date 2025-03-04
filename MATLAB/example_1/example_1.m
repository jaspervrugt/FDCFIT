%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FFFFFFFFFFF DDDDDDDDD  CCCCCCCCCC FFFFFFFFFFF IIIIIIIIIIII TTTTTTTTTT  %
%  FFFFFFFFFFF DDDDDDDDDD CCCCCCCCC  FFFFFFFFFFF  IIIIIIIIII  TTTTTTTTTT  %
%  FF          DD      DD CC         FF               II          TT      %
%  FF          DD      DD CC         FF               II          TT      %
%  FFFFFF      DD      DD CC         FFFFFF           II          TT      %
%  FF          DD      DD CC         FF               II          TT      %
%  FF          DDDDDDDDDD CCCCCCCCC  FF           IIIIIIIIII      TT      %
%  FF          DDDDDDDDD  CCCCCCCCCC FF          IIIIIIIIIIII     TT      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE        1111   %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE           11 11   %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       11  11   %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE              11   %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE          11   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Following models are part of the FDCFIT toolbox
% model       name      Acronym used in paper Sadegh, Vrugt et al., 2015
%   1   'lognormal-2'   LN-2
%   2   'gumbel'        G
%   3   'logistic'      LG
%   4   'logarithmic'   LOG
%   5   'power'         PW
%   6   'quimpo'        Q
%   7   'viola'         V
%   8   'genuchten-2'   VG-2
%   9   'kosugi-2'      K-2
%   10  'lognormal-3'   LN-3
%   11  'pareto'        GP
%   12  'gev'           GEV
%   13  'franchini'     F
%   14  'genuchten-3'   VG-3
%   15  'kosugi-3'      K-3

% clear workspace and figures
clc; clear; close all hidden;      

% Define model and formulation
model = 10;         % Which model do you want to use? Input: [1...15]
FDCPar.form = 'e';  % Fitting approach (fit to exceedance probabilities)
method = 'de';      % Optimization method 
                    % 'LM'    Levenberg Marquardt          (= local method)
                    % 'SP'    Nelder-Mead Simplex          (= local method)
                    % 'CMA'   Covariance Matrix Adaptation (= global methd)
                    % 'DE'    Differential Evolution       (= global methd)

% Define fields of structure options
options.N = 5;              % # trials with optimization method
options.TolX = 1e-2;        % Termination criteria on parameters                                  
options.TolFun = 1e-3;      % Termination criteria objective function (SSE)
options.MaxFunEvals = 1e4;  % Maximum # function evaluations each trial 
options.type = 'day';       % Time scale of flow duration curve
options.print = 'yes';      % Output to screen ('yes' or 'no')

% Example 1: French Broad River - wettest of MOPEX data set
ID_watershed = '08167500';

% Now compute the daily empirical FDC (discharge in column 6)
[e,y,p_0] = calc_FDC(ID_watershed,options.type,6);

% Run FDCFIT toolbox
[x,RMSE] = FDCFIT(FDCPar,model,e,y,method,options);

