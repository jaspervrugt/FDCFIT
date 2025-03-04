function [X,FDCPar,Par_info,options,lsq_options,n_crash,units, ...
    math_str,count] = FDCFIT_setup(FDCPar,model,e,method,options,f_names)
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
% Defines parameter ranges for each model class and formulation           %
%                                                                         %
% SYNOPSIS: [X,FDCPar,Par_info,options,lsq_options,index,units, ...       %
%               math_str,count] = FDCFIT_setup(FDCPar,model,e,method, ... %
%               options,f_names)                                          %
%  where                                                                  %
%   FDCPar    [input] Structure with settings for FDCFIT                  %
%   model     [input] Integer of parametric FDC expression                %
%      = 1          lognormal-2'                                          %
%      = 2          gumbel                                                %
%      = 3          logistic                                              %
%      = 4          logarithmic                                           %
%      = 5          power                                                 %
%      = 6          quimpo                                                %
%      = 7          viola                                                 %
%      = 8          genuchten-2                                           %
%      = 9          kosugi-2                                              %
%      = 10         lognormal-3                                           %
%      = 11         pareto                                                %
%      = 12         gev                                                   %
%      = 13         franchini                                             %
%      = 14         genuchten-3                                           %
%      = 15         kosugi-3                                              %
%      = 16         weibull                                               %
%      = 17         non_central_F                                         %
%      = 18         exponential                                           %
%      = 19         gamma                                                 %
%   e         [input] nx1 vector of empirical exceedance probabilities    %
%   method    [input] Name (string) of parameter estimation method        %
%    = 'lm'         Levenberg-Marquardt method                            %
%    = 'sp'         Nelder-Mead simplex algorithm                         %
%    = 'de'         Differential evolution                                %
%    = 'cma'        Covariance-matrix adaptation                          %
%   options   [input] Optional structure for algorithmic variables        %
%   f_names   [input] Name (string) of suite FDC parametric expressions   %
%   X         [outpt] initial draws of parameter values                   %
%   etc.                                                                  %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Dec. 2014                               %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Initialize variable
n_crash = 0; count = 0;
% open an output file with warnings
fid = fopen('warning_file.txt','a+','n');
% Define model_names
model_acronyms = {'LN','G','LG','LOG','PW','Q','V','VG','K', ...
    'LN','GP','GEV','FS','VG','K','WBL','NCF','EXP','GAM'};
% Define number of parameters of each model
n_pars = [ 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 2 3 2 2 ];
% Now determine how many parameters?
FDCPar.d = n_pars(model);
% User has specified model to be integer values --> change to name model
FDCPar.acronym = char(model_acronyms(model));
% Store model name as lower
FDCPar.model = char(f_names(model));
% Store number of data points
FDCPar.n = numel(e);

% Calculate E_star --> input to many FDC expressions
% Meas_info.E_star = Meas_info.E / (1 - p_0);

% Define the maximum values for the parameters of each model and formulation
switch FDCPar.model
    case    'lognormal-2'
        min = [ -inf 0 ] ; max = [ inf inf ]; units = {'(L/T)','(L/T)'};
        if strcmp(FDCPar.form,'y')
            math_str = ['${\rm Lognormal-2:} \hspace{2mm} ' ...
                'y = \exp\left(a_{\rm LN} - \sqrt{2}b_{\rm LN} ' ...
                '{\rm erfc}^{-1}\bigl(2(1 - e_{\rm n})\bigr)\right)$'];
        else
            math_str = ['${\rm Lognormal-2:} \hspace{2mm} ' ...
                'e_{\rm n} = 1-\frac{1}{2}{\rm erfc} ' ...
                '\left(\frac{a_{\rm LN} - \log(y)}{\sqrt{2} ' ...
                '\: b_{\rm LN}}\right)$'];
        end
    case    'gumbel'
        min = [ -inf 0 ] ; max = [ inf inf ]; units = {'(L/T)','(L/T)'};
        if strcmp(FDCPar.form,'y')
            math_str = ['${\rm Gumbel:} \hspace{2mm} ' ...
                'y = a_{\rm G} - b_{\rm G}\log\left[\log' ...
                '\left(\frac{1}{(1 - e_{\rm n})}\right) \right]$'];
        else
            math_str = ['${\rm Gumbel:} \hspace{2mm} ' ...
                'e_{\rm n} =  1 - \exp\left[-\exp\left(\frac' ...
                '{a_{\rm G} - y}{b_{\rm G}}\right)\right] $'];
        end
    case    'logistic'
        min = [ -inf 0 ] ; max = [ inf inf ]; units = {'(L/T)','(L/T)'};
        if strcmp(FDCPar.form,'y')
            math_str = ['${\rm Logistic:}  \hspace{2mm} ' ...
                'y = a_{\rm LG} - b_{\rm LG}\log\left(' ...
                '\frac{1}{(1-e_{\rm n})} - 1\right)$'];
        else
            math_str = ['${\rm Logistic:} \hspace{2mm} e_{\rm n} = ' ...
                '1 - \left[1 + \exp\left(\frac{a_{\rm LG} - y}' ...
                '{b_{\rm LG}}\right)\right]^{-1}$'];
        end
    case    'logarithmic'
        min = [ -inf -inf ] ; max = [ 0 inf ]; units = {'(L/T)','(L/T)'};
        if strcmp(FDCPar.form,'y')
            math_str = ['${\rm Logarithmic:} \hspace{2mm} ' ...
                'y = b_{\rm LOG} + a_{\rm LOG}\log(e_{\rm n})$'];
        else
            math_str = ['${\rm Logarithmic:} \hspace{2mm} ' ...
                'e_{\rm n} = \exp\left( \frac{y - b_{\rm LOG}}' ...
                '{a_{\rm LOG}}\right) $'];
        end
    case    'power'
        min = [ 0 0 ] ; max = [ inf inf ]; units = {'(-)','(L/T)'};
        if strcmp(FDCPar.form,'y')
            math_str = ['${\rm Power:} \hspace{2mm} ' ...
                'y = b_{\rm PW}e_{\rm n}^{- a_{\rm PW} }$'];
        else
            math_str = ['${\rm Power:} \hspace{2mm} ' ...
                'e_{\rm n} = \left(\frac{1}{b_{\rm PW}}y\right)' ...
                '^{(-1/a_{\rm PW})} $'];
        end
    case    'quimpo'
        min = [ 0 -inf ] ; max = [ inf inf ]; units = {'(L/T)','(-)'};
        if strcmp(FDCPar.form,'y')
            math_str = ['${\rm Quimpo:} \hspace{2mm} ' ...
                'y = a_{\rm Q} \exp\left(-b_{\rm Q}e_{\rm n}\right)$'];
        else
            math_str = ['${\rm Quimpo:} \hspace{2mm} ' ...
                'e_{\rm n} = - \frac{1}{b_{\rm Q}}\log' ...
                '\left(\frac{1}{a_{\rm Q}}y\right) $'];
        end
    case    'viola'
        min = [ 0 0 ] ; max = [ inf inf ]; units = {'(L/T)','(-)'};
        if strcmp(FDCPar.form,'y')
            math_str = ['${\rm Viola:} \hspace{2mm} ' ...
                'y = a_{\rm V} \left(\frac{1-e_{\rm n}}' ...
                '{e_{\rm n}}\right)^{b_{\rm V}}$'];
        else
            math_str = ['${\rm Viola:} \hspace{2mm}' ...
                'e_{\rm n} = \left[\left(\frac{y}{a_{\rm V}}' ...
                '\right)^{(1/b_{\rm V})} + 1\right]^{-1} $'];
        end
    case    'genuchten-2'
        min = [ 0 0 ] ; max = [ inf inf ]; units = {'(T/L)','(-)'};
        if strcmp(FDCPar.form,'y')
            math_str = ['${\rm Genuchten-2:} \hspace{2mm} ' ...
                'y = \frac{1}{a_{\rm VG}} \left[e_{\rm n}' ...
                '^{\bigl(-b_{\rm VG}/(1-b_{\rm VG})\bigr)} - 1\right]' ...
                '^{(1/b_{\rm VG})}$'];
        else
            math_str = ['${\rm Genuchten-2:} \hspace{2mm} ' ...
                'e_{\rm n} = \left[1 + (a_{\rm VG}{y})^{b_{\rm VG}}' ...
                '\right]^{(1/b_{\rm VG}-1)}$'];
        end
    case    'kosugi-2'
        min = [ 0 0 ] ; max = [ inf inf ]; units = {'(L/T)','(-)'};
        if strcmp(FDCPar.form,'y')
            math_str = ['${\rm Kosugi-2:} \hspace{2mm} ' ...
                'y = a_{\rm K}\exp\left(\sqrt{2}b_{\rm K}' ...
                '{\rm erfc}^{-1}(2e_{\rm n})\right)$'];
        else
            math_str = ['${\rm Kosugi-2:} \hspace{2mm} ' ...
                'e_{\rm n} = \frac{1}{2}{\rm erfc}' ...
                '\left[\frac{1}{\sqrt{2} b_{\rm K}}' ...
                '\log\left(\frac{y}{a_{\rm K}}\right)\right]$'];
        end
    case    'lognormal-3'
        min = [ -inf 0 -inf ] ; max = [ inf inf inf ]; 
        units = {'(L/T)','(L/T)','(L/T)'};
        if strcmp(FDCPar.form,'y')
            math_str = ['${\rm Lognormal-3:} \hspace{2mm} ' ...
                'y = c_{\rm LN} + \exp\left[ a_{\rm LN} - \sqrt{2}' ...
                'b_{\rm LN}{\rm erfc}^{-1}\bigl( 2(1-e_{\rm n})' ...
                '\bigr)\right]$'];
        else
            math_str = ['${\rm Lognormal-3:} \hspace{2mm} ' ...
                'e_{\rm n} = \left\{\begin{array}{ll} 1 - \frac{1}{2}' ...
                '{\rm erfc}\left[\frac{a_{\rm LN} - \log(y - ' ...
                'c_{\rm LN})}{\sqrt{2} b_{\rm LN}}\right] ' ...
                '\hspace{0.5cm} {\rm if} \quad y > c_{\rm LN}' ...
                '\\ 1\hphantom{- \frac{1}{2}' ...
                '{\rm erfc}\left[\frac{a_{\rm LN} - \log(y - ' ...
                'c_{\rm LN})}{\sqrt{2} b_{\rm LN}}\right]} ' ...
                '\hspace{0.7cm} {\rm if} \quad y \leq c_{\rm LN}' ...
                '\end{array} \right.$'];
        end
    case    'pareto'
        min = [ -inf 0 -inf ] ; max = [ inf inf inf ]; 
        units = {'(L/T)','(L/T)','(-)'};
        if strcmp(FDCPar.form,'y')
            math_str = ['${\rm Pareto:} \hspace{2mm} y = a_{\rm GP} + ' ...
                '\frac{b_{\rm GP}}{c_{\rm GP} } \left[e_{\rm n}^' ...
                '{-c_{\rm GP}} - 1 \right]$'];
        else
            math_str = ['${\rm Pareto:} \hspace{2mm} ' ...
                'e_{\rm n} = 1 - \left \{\begin{array}{ll} ' ...
                '\left[1+ \frac{c_{\rm GP}(y - a_{\rm GP})}' ...
                '{b_{\rm GP}}\right]^{(-1/c_{\rm GP})} ' ...
                '\hspace{0.16cm} {\rm if} \quad c_{\rm GP} \neq 0 ' ...
                '\\ \exp\left(\frac{a_{\rm GP} - y}{b_{\rm GP}}' ...
                '\right) \hspace{1.7cm} {\rm if} \quad c_{\rm GP} = 0' ...
                '\end{array}\right.$'];
        end
    case    'gev'
        min = [ -inf 0 -inf ] ; max = [ inf inf inf ]; 
        units = {'(L/T)','(L/T)','(-)'};
        if strcmp(FDCPar.form,'y')
            math_str = ['${\rm GEV:} \hspace{2mm} y = a_{\rm GEV} + ' ...
                '\frac{b_{\rm GEV}}{c_{\rm GEV}}\left[\log' ...
                '\left(\frac{1}{(1 - e_{\rm n})}\right)^' ...
                '{-c_{\rm GEV}} - 1\right]$'];
        else
            math_str = ['${\rm GEV:} \hspace{2mm} e_{\rm n} = ' ...
                '\left\{\begin{array}{ll} 1 - \exp\left\{' ...
                '-\left[1+c_{\rm GEV}\left(\frac{y- a_{\rm GEV}}' ...
                '{b_{\rm GEV}}\right)\right]^{(-1/c_{\rm GEV})}' ...
                '\right\} \hspace{0.1cm} {\rm if}\quad c_{\rm GEV}' ...
                ' \neq 0  \\ 1 - \exp \left[ -\exp \left(' ...
                '\frac{ a_{\rm GEV} - y }{ b_{\rm GEV} }\right)' ...
                '\right] \hspace{2.59cm} {\rm if} \quad ' ...
                'c_{\rm GEV} = 0 \end{array} \right. $'];
        end
    case    'franchini'
        min = [ 0 -inf 0 ] ; max = [ inf inf inf ]; 
        units = {'(L/T)','(L/T)','(-)'};
        if strcmp(FDCPar.form,'y')
            math_str = ['${\rm Franchini:} \hspace{2mm} ' ...
                'y = b_{\rm FS} + a_{\rm FS}(1-e_{\rm n})^{c_{\rm FS}}$'];
        else
            math_str = ['${\rm Franchini:} \hspace{2mm} ' ...
                'e_{\rm n} = 1 - \left(\frac{ y - b_{\rm FS}}' ...
                '{a_{\rm FS}}\right)^{(1/c_{\rm FS})}$'];
        end
    case    'genuchten-3'
        min = [ 0 1 0 ] ; max = [ inf inf 1 ]; 
        units = {'(T/L)','(-)','(-)'};
        if strcmp(FDCPar.form,'y')
            math_str = ['${\rm Genuchten-3:} \hspace{2mm} ' ...
                'y = \frac{1}{a_{\rm VG}}\left[e_{\rm n}^' ...
                '{(-1/c_{\rm VG})} - 1\right]^{(1/b_{\rm VG})}$'];
        else
            math_str = ['${\rm Genuchten-3:} \hspace{2mm} ' ...
                'e_{\rm n} = \left[ 1 + (a_{\rm VG}{y})^{b_{\rm VG}}' ...
                '\right]^{-c_{\rm VG}}$'];
        end
    case    'kosugi-3'
        min = [ 0 0 -inf ] ; max = [ inf inf inf ]; units = {'(L/T)', ...
            '(-)','(L/T)'};
        if strcmp(FDCPar.form,'y')
            math_str = ['${\rm Kosugi-3:} \hspace{2mm} ' ...
                'y = c_{\rm K} + (a_{\rm K} - c_{\rm K})\exp\left[' ...
                '\sqrt{2}b_{\rm K}\mbox{erfc}^{-1}\left(2e_{\rm n}' ...
                '\right)\right]$'];
        else
            math_str = ['${\rm Kosugi-3:} \hspace{2mm} ' ...
                'e_{\rm n} = \left\{\begin{array}{ll}' ...
                '\frac{1}{2} \: {\rm erfc}\left[\frac{1}{\sqrt{2} ' ...
                '\: b_{\rm K}}\log\left(\frac{y - c_{\rm K} }' ...
                '{ a_{\rm K} - c_{\rm K} }\right)\right] ' ...
                '\hspace{0.5cm} {\rm if}\quad y > c_{\rm K} \\ 1 ' ...
                '\hspace{4.15cm} {\rm if} \quad y \leq c_{\rm  K} ' ...
                '\end{array} \right. $'];
        end
    case    'weibull'
        min = [ 0 0 ] ; max = [ inf inf ]; math_str = [];
    case    'non_central_F'
        min = [ 0 0 0 ] ; max = [ inf inf inf ]; math_str = [];
end

% Define default ranges
Par_info.min = min; Par_info.max = max;

% Check settings
if ~isfield(options,'TolX')
    % Default setting of number of trials
    options.TolX = 1e-2;
end

% Check settings
if ~isfield(options,'TolFun')
    % Default setting of number of trials
    options.TolFun = 1e-3;
end

% Check SP/LM settings
if ~isfield(options,'MaxFunEvals')
    % Default setting of number of trials
    options.MaxFunEvals = 1e4;
end

% Check SP/LM settings
if any(strcmpi(method,{'sp','lm'}))
    if ~isfield(options,'N')
        % Default setting of number of trials
        options.N = 5;
    end
end
% Verify DE/CMAES settings
if any(strcmpi(method,{'de','cma'}))
    if ~isfield(options,'P')
        % Default setting of population size
        options.P = 25;
    end
end

% Verify DE settings
if strcmpi(method,'de')
    if ~isfield(options,'CR')
        % Default setting of crossover used by DE
        options.CR = 0.8;
    end
end

if any(strcmpi(method,{'de','cmaes'}))
    lsq_options = [];
else
    % Define optim options
    lsq_options = optimset('Display','off','TolX', ...
        options.TolX,'TolFun',options.TolFun, ...
        'MaxFunEvals',options.MaxFunEvals);
end

% Now check this
if ~any(strcmpi(method,{'de','cma'}))
    % Remove field P - not needed
    if isfield(options,'P')
        options = rmfield(options,'P');
    end
end

% Now check this
if ~any(strcmpi(method,{'lm','sp'}))
    % Remove field P - not needed
    if isfield(options,'N')
        options = rmfield(options,'N');
    end
end

if ~strcmp(method,'de')
    if isfield(options,'CR')
        options = rmfield(options,'CR');
    end
    if isfield(options,'F')
        options = rmfield(options,'F');
    end
end

if ~isfield(options,'print')
    options.print = 'yes';
end

if ~isfield(options,'type')
    options.type = 'daily';
end

% Now generate prior ranges and draw samples
max(isinf(max)) = 2; min(isinf(min)) = -2;

% Create samples for multi-start method or for DE
switch method
    case 'de'   % Draw initial population DE
        X = LH_sampling(min,max,options.P); 
    case 'lm'   % Draw initial population LM
        X = LH_sampling(min,max,options.N); 
    case 'sp'   % Draw multistart points SP
        X = LH_sampling(min,max,options.N); 
    otherwise   % Empty for CMAES (creates own)
        X = [];                         
end

% Now print to screen all the settings
disp('---------------- Summary of fields FDCPar -----------------');
disp(FDCPar)
disp('-----------------------------------------------------------');

% Now print this only if DE is used
disp('--------------- Summary of fields options -----------------');
disp(options)
disp('-----------------------------------------------------------');

% Now close warning_file.txt file
fclose(fid);    

end
