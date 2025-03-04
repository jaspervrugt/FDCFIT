% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  FFFFFFFFFFF DDDDDDDDD  CCCCCCCCCC FFFFFFFFFFF IIIIIIIIIIII TTTTTTTTTT  %
%  FFFFFFFFFFF DDDDDDDDDD CCCCCCCCC  FFFFFFFFFFF  IIIIIIIIII  TTTTTTTTTT  %
%  FF          DD      DD CC         FF               II          TT      %
%  FF          DD      DD CC         FF               II          TT      %
%  FFFFFF      DD      DD CC         FFFFFF           II          TT      %
%  FF          DD      DD CC         FF               II          TT      %
%  FF          DDDDDDDDDD CCCCCCCCC  FF           IIIIIIIIII      TT      %
%  FF          DDDDDDDDD  CCCCCCCCCC FF          IIIIIIIIIIII     TT      %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
% This code determines the fitting coefficients of a large suite of       %
% models of the flow duration curve. This includes a new set of           %
% parametric expressions based on the water retention functions of van    %
% Genuchten (1980) and Kosugi (1996).                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% SYNOPSIS: [x,RMSE] = FDCFIT(FDCPar,model,e,y,method)                    %
%           [x,RMSE] = FDCFIT(FDCPar,model,e,y,method,options)            %
%  where                                                                  %
%   FDCPar    [input] Structure with settings for FDCFIT                  %
%    .d           # unknown parameters of parametric FDC expression       %
%    .form        Forward or inverse formulation of the FDC               %
%      = 'e'      Forward expression: e_s = f(y)                          %
%      = 'y'      Backward/inverse expression: y_s = f(e)                 %
%    .n           # discharge measurements [= numel(y)]  return by FDCFIT %
%    .model       Full name of FDC expression            return by FDCFIT %
%    .acronym     Acroynym FDC expression                return by FDCFIT %
%   model     [input] Integer of parametric FDC expression                %
%      = 1          Full name: lognormal-2         Acronym: LN-2          %
%      = 2          Full name: gumbel              Acronym: G             %
%      = 3          Full name: logistic            Acronym: LG            %
%      = 4          Full name: logarithmic         Acronym: LOG           %
%      = 5          Full name: power               Acronym: PW            %
%      = 6          Full name: quimpo              Acronym: Q             %
%      = 7          Full name: viola               Acronym: V             %
%      = 8          Full name: genuchten-2         Acronym: VG-2          %
%      = 9          Full name: kosugi-2            Acronym: K-2           %
%      = 10         Full name: lognormal-3         Acronym: LN-3          %
%      = 11         Full name: pareto              Acronym: GP            %
%      = 12         Full name: gev                 Acronym: GEV           %
%      = 13         Full name: franchini           Acronym: F             %
%      = 14         Full name: genuchten-3         Acronym: VG-3          %
%      = 15         Full name: kosugi-3            Acronym: K-3           %
%   e         [input] nx1 vector of empirical exceedance probabilities    %
%   y         [input] nx1 vector of empirical (= measured) discharges     %
%   method    [input] Name (string) of parameter estimation method        %
%    = 'lm'         Levenberg-Marquardt method                            %
%    = 'sp'         Nelder-Mead simplex algorithm                         %
%    = 'de'         Differential evolution                                %
%    = 'cma'        Covariance-matrix adaptation                          %
%   options   [input] Optional structure for algorithmic variables        %
%    .N             # trials 'lm' and 'sp' methods        DEF: 5          %
%    .P             Population size for 'de' and 'cma'    DEF: 25         %
%    .CR            Crossover value for 'de' algorithm    DEF: 0.8        %
%    .TolX          Error tolerance on parameter values   DEF: 1e-3       %
%    .TolFun        Error tolerance on objective function DEF: 1e-4       %
%    .MaxFunEvals   Maximum # function evaluations        DEF: 1e4        %
%    .type          Type of FDC derived from daily discharge data, y      %
%     = 'day'       Daily flow duration curve                             %
%     = 'week'      Weekly flow duration curve                            %
%     = 'month'     Monthly flow duration curve                           %
%     = 'annual'    Annual flow duration curve                            %
%    .print         Output writing screen (tables/figs)   DEF: 'yes'      %
%     = 'no'        No output writing to the screen                       %
%     = 'yes'       Output writing to the screen          DEFault         %
%   x         [outpt] Maximum likelihood parameters FDC expressions       %
%   RMSE      [outpt] Root Mean Square Error of FDC expression            %
%                                                                         %
% LITERATURE:                                                             %
%  Vrugt, J.A. (2017), FDCFIT: A MATLAB toolbox of parametric             %
%      expressions of the Flow duration curve, Manual, pp. 1-37.          %
%  Sadegh, M., J.A. Vrugt, X. Cu, and H.V. Gupta, (2016), The soil water  %
%      characteristic as new class of parametric expressions of the flow  %
%      duration curve, Journal of Hydrology, 535, pp. 438-456,            %
%          https://doi.org/10.1016/j.jhydrol.2016.01.027                  %
%  Vrugt, J.A., and M. Sadegh (2013), Toward diagnostic model             %
%      calibration and evaluation: Approximate Bayesian computation,      %
%      Water Resources Research, 49, 4335–4345,                           %
%          https://doi.org/10.1002/wrcr.20354                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% (c) Written by Jasper A. Vrugt, Dec. 2014                               %
% University of California Irvine                                         %
%  Version 1    December, 2014                                            %
%  Version 1.1  Feb, 2016                                                 %
%  Version 1.2  May, 2017                                                 %
%  Version 2    July, 2024                                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add main FDCFIT directory and underlying postprocessing directory
addpath(pwd,[pwd,'\miscellaneous']);
% Now go to example 1
cd example_1
% And you can now execute example_1 using (uncomment)
% example_1
