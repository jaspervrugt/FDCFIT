function [FDCPar,options,f_names] = FDCFIT_check(FDCPar,model,e,y, ...
    method,options)
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
% Checks the settings defined by the user                                 %
%                                                                         %
% SYNOPSIS: [FDCPar,options,f_names] = FDCFIT_check(FDCPar,model,e,y, ... %
%               method,options)                                           %
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
%   y         [input] nx1 vector of empirical (= measured) discharges     %
%   method    [input] Name (string) of parameter estimation method        %
%    = 'lm'         Levenberg-Marquardt method                            %
%    = 'sp'         Nelder-Mead simplex algorithm                         %
%    = 'de'         Differential evolution                                %
%    = 'cma'        Covariance-matrix adaptation                          %
%   options   [input] Optional structure for algorithmic variables        %
%   FDCPar    [outpt] Structure with settings for FDCFIT [verified]       %
%   options   [outpt] Optional structure algorithmic variables [verified] %
%   f_names   [outpt] Name (string) of suite FDC parametric expressions   %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Dec. 2014                               %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


% % Derive current time and set deadline
% deadline = datenum('28-Feb-2025');
% 
% % Now check whether this is a trial version or not
% if ( deadline - now ) < 0
%     % ERROR -- trial version ended
%     error('FDCFIT ERROR: Trial version of FDCFIT V1.1 has ended');
% end

% open an output file with warnings
fid = fopen('warning_file.txt','w+');
fprintf(fid,'-------------- FDCFIT warning file --------------\n');

% Check FDCFIT input arguments
if isempty(FDCPar)
    error(['FDCFIT ERROR: input argument FDCPar should not be empty ' ...
        'but a structure with fields']);
elseif ~isstruct(FDCPar)
    error(['FDCFIT ERROR: input argument FDCPar should be a structure ' ...
        'with fields']);
end
% Check FDCFIT input arguments
if isempty(model)
    error(['FDCFIT ERROR: input argument model should not be empty ' ...
        'but an integer of [1,2,...,15]']);
elseif ~isnumeric(model)
    error(['FDCFIT ERROR: input argument model should not be empty ' ...
        'but an integer of [1,2,...,15]']);
end
% Check FDCFIT input arguments
if isempty(e)
    error(['FDCFIT ERROR: input argument e should not be empty but ' ...
        'a nx1 vector of empirical exceedance probabilities']);
elseif ~isnumeric(e)
    error(['FDCFIT ERROR: input argument e should not be empty but ' ...
        'a nx1 vector of empirical exceedance probabilities']);
end
% Check FDCFIT input arguments
if isempty(y)
    error(['FDCFIT ERROR: input argument y should not be empty but ' ...
        'a nx1 vector of measured discharges']);
elseif ~isnumeric(y)
    error(['FDCFIT ERROR: input argument y should not be empty but ' ...
        'a nx1 vector of measured discharges']);
end
% Check FDCFIT input arguments
if isempty(options)
    error(['FDCFIT ERROR: input argument options should not be ' ...
        'empty but a structure with fields']);
elseif ~isstruct(options)
    error(['FDCFIT ERROR: input argument options should be a ' ...
        'structure with fields']);
end

% Define model_names
f_names = {'lognormal-2'
    'gumbel'
    'logistic'
    'logarithmic'
    'power'
    'quimpo'
    'viola'
    'genuchten-2'
    'kosugi-2'
    'lognormal-3'
    'pareto'
    'gev'
    'franchini'
    'genuchten-3'
    'kosugi-3'
    'weibull'
    'non_central_F'
    'exponential'
    'gamma'};

%<><><><><><><><><><><><><><><><> FDCPar <><><><><><><><><><><><><><><><><>
if ~isfield(FDCPar,'form')
    % Now print to screen
    error(['FDCFIT ERROR: The field ''form'' of structure FDCPar ' ...
        'is not defined: Set its content to ''e'' or ''y'' respectively']);
elseif isempty(FDCPar.form)
    % expression should be natural numbers
    error(['FDCFIT ERROR: Field ''formulation'' of structure FDCPar ' ...
        'should not be empty but set to ''e'' or ''y'' respectively']);
elseif ~char(FDCPar.form)
    % expression should be natural numbers
    error(['FDCFIT ERROR: Field ''form'' of structure FDCPar should ' ...
        'list a string enclosed between quotes: ''e'' or ''y''' ...
        ' respectively']);
else
    % First make it lower case letters
    FDCPar.form = lower(FDCPar.form);
    % Now check
    if ~any(strcmp(FDCPar.form,{'e','y'}))
        % Now print to screen
        error(['FDCFIT ERROR: Unknown option used for field ' ...
            '''form'' of structure FDCPar - should be ''e'' or ''y''']);
    end
end

%<><><><><><><><><><><><><><><><> model <><><><><><><><><><><><><><><><><><
% Verify input settings
if isempty(model)
    % expression should be natural numbers
    error(['FDCFIT ERROR: Input argument ''model'' should not be ' ...
        'empty but list an integer from [1,2,...,15]']);
elseif ~isnumeric(model)
    % expression should be natural numbers
    error(['FDCFIT ERROR: Input argument ''model'' should not be ' ...
        'a string but an integer selected from [1,2,...,15]']);
end
if isnumeric(model)
    % Check whether rounded
    if mod(model,round(model))
        % expression should be natural numbers
        error(['FDCFIT ERROR: Input argument ''model'' should be ' ...
            'an integer selected from [1,2,...,15]']);
    end
    % Now check again
    if (model < 1) || (model > 15)
        % expression should be natural numbers
        error(['FDCFIT ERROR: Inpurt argument ''model'' should be ' ...
            'an integer selected from [1,2,...,15]']);
    end
else                            % user has specified string
    % Print available models
    for i = 1 : 15
        fprintf(' [ %1d ] %s \n',i,char(f_names(i)));
    end
    % Now print to screen
    error(['FDCFIT ERROR: Input argument ''model'' should be ' ...
        'an integer from [1,2,...,15] (see listed models above)']);
end

%<><><><><><><><><><><><><><><><> e and y <><><><><><><><><><><><><><><><><
% Compute size of e and y
a = size(e); b = size(y);
if sum(a > 1) > 1
    % Now print to screen
    error(['FDCFIT ERROR: Input argument ''e'' should not be a matrix ' ...
        'but a nx1 vector']);
end
if sum(b > 1) > 1
    % Now print to screen
    error(['FDCFIT ERROR: Input argument ''y'' should not be a matrix ' ...
        'but a nx1 vector']);
end
% Verify input settings
if a(1) ~= b(1)
    % Now print to screen
    error(['FDCFIT ERROR: Input arguments ''y'' and ''e'' do not ' ...
        'match in size']);
end
if a(1) < 3
    % Now print to screen
    error(['FDCFIT ERROR: Insufficient number of discharge/exceedance ' ...
        'probability observations']);
end
if a(1) < 25
    % Now print to screen
    evalstr = char(['FDCFIT WARNING: Rather small number of ' ...
        'observations of discharge/exceedance probability\n']);
    % Now print warning to screen and to file
    fprintf(evalstr); fprintf(fid,evalstr);
end

%<><><><><><><><><><><><><><><><> method <><><><><><><><><><><><><><><><><>
if isempty(method)
    % expression should be natural numbers
    error(['FDCFIT ERROR: Input argument ''method'' should not be ' ...
        'empty but list ''SP'' or ''LM'' or ''DE'' or ''CMA'' ' ...
        'respectively']);
elseif ~char(method)
    % expression should be natural numbers
    error(['FDCFIT ERROR: Input argument ''method'' should list a ' ...
        'string enclosed between quotes: ''SP'' or ''LM'' or ''DE'' ' ...
        'or ''CMA'' respectively']);
else
    % Now check
    if ~any(strcmpi(method,{'lm','sp','de','cma'}))
        % ERROR -- unknown method
        error(['FDCFIT ERROR: Unknown optimization method ' ...
            '-> Set method = ''LM'' or ''SP'' or ''DE'' or ''CMA'' ']);
    end
end
if any(strcmp(method,{'lm','sp'}))
    % Check SP/LM settings
    if ~isfield(options,'N')
        % Print warning
        evalstr = char(strcat(['FDCFIT WARNING: Field ''N'' (number ' ...
            'trials with SP/LM) of structure options not defined ' ...
            'by user --> options.N = 5 ( default value )\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Remove field
        options = rmfield(options,'N');
    elseif isempty(options.N)
        % Print warning
        evalstr = char(strcat(['FDCFIT WARNING: Field ''N'' (number ' ...
            'trials with SP/LM) of structure options is left empty ' ...
            'by user --> options.N = 5 ( default value )\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Remove field
        options = rmfield(options,'N');
    elseif ~isnumeric(options.N)
        % Print warning
        evalstr = char(strcat(['FDCFIT WARNING: Field ''N'' (number ' ...
            'trials with SP/LM) of structure options should be an ' ...
            'integer --> options.N = 5 ( default value )\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Remove field
        options = rmfield(options,'N');
    else
        if options.N < 0
            % Now print to screen
            evalstr = char(strcat(['FDCFIT WARNING: The field ''N'' ' ...
                'of structure options cannot be negative (default ' ...
                'is options.N = 5)\n']));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);            
            % Remove field
            options = rmfield(options,'N');            
        end
        if options.N < 5
            evalstr = char(strcat('FDCFIT WARNING: The value of',{' '}, ...
                num2str(options.N),{' '},['different trials with LM ' ...
                'or SP is rather small (default is options.N = 5)\n']));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
        end
        if options.N > 25
            evalstr = char(strcat('FDCFIT WARNING: The value of', ...
                {' '},num2str(options.N),{' '},['different trials ' ...
                'with LM or SP is rather large ' ...
                '(default is options.N = 5)\n']));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
        end
        if mod(options.N,round(options.N))
            % Print warning
            evalstr = char(strcat(['FDCFIT WARNING: Field ''N'' ' ...
                '(number trials with SP/LM) of structure options ' ...
                'should be an integer ' ...
                '--> options.N = 5 (default value)\n']));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
            % Remove field
            options = rmfield(options,'N');
        end
    end
end
%<><><><><><><><><><><><><><><><> options <><><><><><><><><><><><><><><><><
% Check whether user has supplied right info about structure options
if ~isfield(options,'TolX')
    % Print warning
    evalstr = char(strcat(['FDCFIT WARNING: Field ''TolX'' of ' ...
        'structure options not defined by user ' ...
        '--> options.TolX = 1e-2 ( default value )\n']));
    % Now print warning to screen and to file
    fprintf(evalstr); fprintf(fid,evalstr);
elseif isfield(options,'TolX')
    if isempty(options.TolX)
        % Print warning
        evalstr = char(strcat(['FDCFIT WARNING: Field ''TolX'' of ' ...
            'structure options is left empty by user ' ...
            '--> default settings options.TolX = 1e-2\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Remove field
        options = rmfield(options,'TolX');
    elseif ~isnumeric(options.TolX)
        % Now print to screen
        evalstr = char(strcat(['FDCFIT WARNING: The field ''TolX'' of ' ...
            'structure optim should be a numerical value ' ...
            '--> resort to default setting options.TolX = 1e-2\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Remove field
        options = rmfield(options,'TolX');
    else
        if options.TolX < 0
            % Now print to screen
            evalstr = char(strcat(['FDCFIT WARNING: The field ' ...
                '''TolX'' of structure options cannot be negative ' ...
                '(default is options.TolX = 1e-2)\n']));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
            % Remove field
            options = rmfield(options,'TolX');            
        end
        if options.TolX > 0.1
            % Now print to screen
            evalstr = char(strcat(['FDCFIT WARNING: The value ' ...
                'of'],{' '},num2str(options.TolX),{' '},['for ' ...
                'options.TolX is rather large ' ...
                '(default is options.TolX = 1e-2)\n']));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
        end
    end
end
% Check settings
if ~isfield(options,'TolFun')
    % Print warning
    evalstr = char(strcat(['FDCFIT WARNING: Field ''TolFun'' of ' ...
        'structure options not defined by user ' ...
        '--> options.TolFun = 1e-4 ( default value )\n']));
    % Now print warning to screen and to file
    fprintf(evalstr); fprintf(fid,evalstr);
elseif isfield(options,'TolFun')
    if isempty(options.TolFun)
        % Print warning
        evalstr = char(strcat(['FDCFIT WARNING: Field ''TolFun'' of ' ...
            'structure options is left empty by user ' ...
            '--> default settings options.TolFun = 1e-3\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Remove field
        options = rmfield(optim,'TolFun');
    elseif ~isnumeric(options.TolFun)
        % Now print to screen
        evalstr = char(strcat(['FDCFIT WARNING: The field ''TolFun'' ' ...
            'of structure options should be a numerical value ' ...
            '--> resort to default setting options.TolFun = 1e-3\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Remove field
        options = rmfield(options,'TolFun');
    else
        if options.TolFun < 0
            % Now print to screen
            evalstr = char(strcat(['FDCFIT WARNING: The field ' ...
                '''TolFun'' of structure options cannot be ' ...
                'negative (default is options.TolFun = 1e-3)\n']));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
            % Remove field
            options = rmfield(options,'TolFun');
        end
        if options.TolFun > 0.1
            % Now print to screen
            evalstr = char(strcat(['FDCFIT WARNING: The value ' ...
                'of'],{' '},num2str(options.TolFun),{' '}, ...
                ['for options.TolFun is rather large ' ...
                '(default is options.TolFun = 1e-3)\n']));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
        end
    end
end
% Check SP/LM settings
if ~isfield(options,'MaxFunEvals')
    % Print warning
    evalstr = char(strcat(['FDCFIT WARNING: Field ''MaxFunEvals'' ' ...
        'of structure options not defined by user ' ...
        '--> options.MaxFunEvals = 1e4 ( default value )\n']));
    % Now print warning to screen and to file
    fprintf(evalstr); fprintf(fid,evalstr);
elseif isfield(options,'MaxFunEvals')
    if isempty(options.MaxFunEvals)
        % Print warning
        evalstr = char(strcat(['FDCFIT WARNING: Field ''MaxFunEvals'' ' ...
            'of structure options is left empty by user ' ...
            '--> default settings options.MaxFunEvals = 1e4\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Remove field
        options = rmfield(options,'MaxFunEvals');
    elseif ~isnumeric(options.MaxFunEvals)
        % Now print to screen
        evalstr = char(strcat(['FDCFIT WARNING: The field ' ...
            '''MaxFunEvals'' of structure options should be a ' ...
            'numerical value ' ...
            '--> resort to default setting options.MaxFunEvals = 1e4\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Remove field
        options= rmfield(options,'MaxFunEvals');
    else
        if options.MaxFunEvals < 0
            % Now print to screen
            evalstr = char(strcat(['FDCFIT WARNING: The field ' ...
                '''MaxFunEvals'' of structure options cannot be ' ...
                'negative (default is options.MaxFunEvals = 1e4)\n']));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
            % Remove field
            options = rmfield(options,'MaxFunEvals');
        end
        if options.MaxFunEvals < 51
            % Now print to screen
            error(['FDCFIT ERROR: The field ''MaxFunEvals'' of ' ...
                'structure options is set too small to provide ' ...
                'reliable results ' ...
                '(default is options.MaxFunEvals = 1e4)']);
        end
        if options.MaxFunEvals > 1e5
            % Now print to screen
            evalstr = char(strcat(['FDCFIT WARNING: The value ' ...
                'of'],{' '},num2str(options.MaxFunEvals),{' '}, ...
                ['for optim.MaxFunEvals is rather large ' ...
                '(default is options.MaxFunEvals = 1e4)\n']));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
        end
    end
end
% Check whether user has supplied right info about FDCPar.options_user
if isfield(options,'type')
    if isempty(options.type)
        % Now print to screen
        evalstr = char(strcat(['FDCFIT WARNING: The field ''type'' of ' ...
            'structure options should not be empty but list a ' ...
            'string (one of ''day''/''week''/''month''' ...
            '/''annual'' )\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    elseif ~ischar(options.type)
        % Now print to screen
        evalstr = char(strcat(['FDCFIT WARNING: The field ''type'' of ' ...
            'structure options should be a string (one of ''day''/' ...
            '''week''/''month''/''annual'' )\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    elseif ~any(strcmp(options.type,{'day','week','month','annual'}))
        % Now print warning to screen
        evalstr = char(strcat(['FDCFIT WARNING: The field ''type'' of ' ...
            'structure options should be a string (one of ''day''/' ...
            '''week''/''month''/''annual'' )']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
end
% Check other fields of options
if isfield(options,'print')
    if isempty(options.print)
        % Now print to screen
        evalstr = char(strcat(['FDCFIT WARNING: The field ''print'' ' ...
            'of structure options should not be empty but list a ' ...
            'string with either ''yes'' or ''no'' ' ...
            '(resort to default ''yes'')\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Remove field
        options = rmfield(options,'print');        
    elseif ~ischar(options.print)
        % Now print to screen
        evalstr = char(strcat(['FDCFIT WARNING: The field ''print'' ' ...
            'of structure options should be a string with ' ...
            'either ''yes'' or ''no'' (resort to default ''yes'')\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
       % Remove field
        options = rmfield(options,'print');           
    elseif ~any(strcmp(options.print,{'no','yes'}))
        % Now print warning to screen
        evalstr = char(strcat(['FDCFIT WARNING: The field ' ...
            '''print'' of structure options should be a string with ' ...
            'either ''yes'' or ''no'' (resort to default ''yes'')\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
       % Remove field
        options = rmfield(options,'print');           
    end
end

% Verify DE/CMAES settings
if any(strcmpi(method,{'de','cma'}))
    if ~isfield(options,'P')
        % Print warning
        evalstr = char(strcat(['FDCFIT WARNING: Field ''P'' ' ...
            '(population size) of structure options not defined ' ...
            'by user --> options.P = 25 ( default value )\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    elseif isempty(options.P)
        % Print warning
        evalstr = char(strcat(['FDCFIT WARNING: Field ''P'' of ' ...
            'structure options is left empty by user ' ...
            '--> default settings options.P = 25\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Remove field
        options = rmfield(options,'P');
    elseif ~isnumeric(options.P)
        % Now print to screen
        evalstr = char(strcat(['FDCFIT WARNING: The field ''P'' ' ...
            'of structure options should be a numerical value ' ...
            '--> resort to default setting options.P = 25\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Remove field
        options = rmfield(options,'P');
    else
        if options.P < 0
            % Now print to screen
            evalstr = char(strcat(['FDCFIT WARNING: The field ' ...
                '''P'' of structure options cannot be negative: ' ...
                'Set to at least ''P = 25'' individuals for ' ...
                'global search with'],{' '},upper(options.method),'\n'));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
            % Remove field
            options.P = rmfield(options,'P');
        end
        if options.P < 10
            % Now print to screen
            evalstr = char(strcat(['FDCFIT WARNING: The field ''P'' ' ...
                'of structure options is set rather small: ' ...
                'At least ''P = 25'' individuals for global ' ...
                'search with'],{' '},upper(options.method),'\n'));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
        end
        if options.P > 50
            % Now print to screen
            evalstr = char(strcat(['FDCFIT WARNING: The field ''P'' ' ...
                'of structure options is set rather large: A ' ...
                'population of ''P = 50'' should be sufficient ' ...
                'for estimating a few parameters with'], ...
                {' '},upper(options.method),'\n'));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
        end
    end
end
% Verify DE settings
if strcmpi(method,'de')
    if ~isfield(options,'CR')
        % Print warning
        evalstr = char(strcat(['FDCFIT WARNING: Field ''CR'' ' ...
            '(crossover value) of structure DE not defined by user ' ...
            '--> options.CR = 0.9 ( default value )\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    elseif isempty(options.CR)
        % Print warning
        evalstr = char(strcat(['FDCFIT WARNING: Field ''CR'' of ' ...
            'structure options is left empty by user ' ...
            '--> default settings options.CR = 0.9\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Remove field
        options = rmfield(options,'CR');
    elseif ~isnumeric(options.CR)
        % Now print to screen
        evalstr = char(strcat(['FDCFIT WARNING: The field ''CR'' of ' ...
            'structure options should be a numerical value ' ...
            '--> resort to default setting options.CR = 0.9\n']));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Remove field
        options = rmfield(options,'CR');
    else
        if options.CR < 0
            % Now print to screen
            evalstr = char(strcat(['FDCFIT WARNING: The field ''CR'' ' ...
                'of structure options cannot be negative: Crossover ' ...
                'value, ''CR'' should be between zero and one\n']));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
            % Remove field
            options = rmfield(options,'CR');
        end
        if options.CR > 1
            % Now print to screen
            error(['FDCFIT ERROR: The field ''CR'' of structure ' ...
                'options cannot be larger than one: Crossover ' ...
                'value, ''CR'' should be between zero and one']);
        end
        if options.CR < 0.4
            % Now print to screen
            error(['FDCFIT ERROR: The field ''CR'' of structure ' ...
                'options is set rather small: ' ...
                'I would recommend to use at least ''CR = 0.6''']);
        end
        if options.CR > 0.9
            % Now print to screen
            error(['FDCFIT ERROR: The field ''CR'' of structure ' ...
                'options is set rather large: A value of ''CR = 0.9'' ' ...
                'would seem better for 2 or 3 parameters']);
        end
    end
end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% Now close warning_file.txt file
fclose(fid);

end
