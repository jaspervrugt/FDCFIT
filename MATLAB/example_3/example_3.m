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
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      333333   %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE              33   %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE          333   %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE              33   %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE      333333   %
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
model = 15;         % Which model do you want to use? Input: [1...15]
FDCPar.form = 'e';  % Fitting approach (fit to exceedance probabilities)
method = 'LM';      % Optimization method 
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

% Example 3: MOPEX data set: watershed_no = [1...10] 
% (unzip file "MOPEX_remaining_watersheds.rar" for all MOPEX data sets)
watershed_no = 10;

% Load the data of basin j
evalstr = strcat('data = load(''data_', ...
    num2str(watershed_no),'.txt'');'); eval(evalstr);
ii = (data==0); 
p_0 = sum(ii)/numel(data);  % Probability of zero flows?
Y = data(data>0);           % Remaining data
Y = Y(~isnan(Y));           % Remove nan due to initialization
N = size(Y,1);              % How many values of Y?
y_s = sort(Y);              % Sort discharge data in ascending order
e_n = 1 - ((1:N)' -0.5)./N; % Exceedance probabilities (based on N only)

% Run FDCFIT toolbox
[x,RMSE] = FDCFIT(FDCPar,model,e_n,y_s,method,options);
