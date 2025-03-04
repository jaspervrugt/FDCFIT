function varargout = FDCFIT_functions(x,FDCPar,e,y,returnarg)
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
% This function evaluates the exceedance probability, e_s, or streamflow, %
% y_s, using a suite of different parametric expressions of the flow      %
% duration curve. Return argument is set by the user by the last input    %
% argument, 1: L2 objective function, 2: residual vector, e-e_s or y-y_s, %
% and 3: L1 objective function                                            % 
%                                                                         %
% SYNOPSIS: varargout = FDCFIT_functions(x,FDCPar,e,y,returnarg)          %
%  where                                                                  %
%   x         [input] 1xd vector parameter values parametric expression   %
%   FDCPar    [input] Structure with settings for FDCFIT                  %
%   e         [input] nx1 vector of empirical exceedance probabilities    %
%   y         [input] nx1 vector of empirical (= measured) discharges     %
%   returnarg [input] Integer of desired return argument                  %
%    = 1            L2 objective function                                 %
%    = 2            Residual vector of exceedance or discharge values     %
%    = 3            L1 objective function                                 %
%   varargout [outpt] L2/L1 objective function or residual vector         %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Dec. 2014                               %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% % % If any parameter < 0 --> return bad SSR
% % if sum(x < 0) > 0
% %     switch returnarg
% %         case 0  % return exceedance probability (LM method)
% %             varargout = 0;
% %         case 1  % return SSR (Simplex method)
% %             varargout = 1e50;
% %         case 2  % return residual (LSQNONLIN)
% %             varargout = 1e10*ones(size(e));      
% %     end
% %     % Now return back to header (bad parameters lead to bad OF values)
% %     return
% % end
a = x(1); b = x(2);

% Which model to use
switch char(FDCPar.model)

    %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                       2 PARAMETER EXPRESSIONS                       %
    %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
    case 'lognormal-2' % lognormal - 2 parameters
        a = x(1); b = x(2);
        switch strcmp(FDCPar.form,'e')
            case 1
                % Calculate log-transformed variable and cdf
                % z = (log(y) - a)./b; cdf = 0.5 * erfc(-z./sqrt(2));
                % Now calculate exceedance probability
                % e_s = e_0 * (1 - cdf); 
                % Now calculate e_s directly
                e_s = 1 - .5 * erfc(-(log(y)-a)/(sqrt(2)*b));            
            otherwise
                % Now write Y as a function of e
                % 1 - e = 0.5 * erfc(-z./sqrt(2));
                % 2*(1-e) = erfc(-z./sqrt(2));
                % erfcinv(2*(1-e)) = -z./sqrt(2);
                % z = -sqrt(2)*erfcinv(2*(1-e));
                % ( log(y_s) - a)./b = -sqrt(2)*erfcinv(2*(1-e));
                % log(y_s) = a - b*sqrt(2)*erfcinv(2*(1-e));
                % y_s = exp(a - b*sqrt(2)*erfcinv(2*(1-e)));
                y_s = exp(a - sqrt(2)*b*erfcinv(2*(1-e)));
        end
    case 'gumbel' % Gumbel - 2 parameters
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % FDC function
                e_s = 1 - exp(-exp(-(y-a)/b));
            otherwise   % Now write Y as function of e
                y_s = a - b * log(-log(1-e));
        end
    case 'logistic' % logistic dist - 2 parameters: http://goo.gl/YRBASP
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % FDC function - after Castellarin 2004 WRR:
                e_s = 1 - 1./(1 + exp(-(y-a)/b));
            otherwise   % Update Y
                y_s = a - b * log(1./(1 - e) - 1);
        end
    case 'logarithmic' % logarithmic - 2 parameters: http://goo.gl/FhQb5o
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % After Sauquet and Catalogne, HESS, 2011
                e_s = exp((y - b)/a);
            otherwise   % Now write Y as a function of e
                y_s = b + a * log(e);
        end
    case 'power' % power - 2 parameters: http://goo.gl/0tfEJw
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % FDC function - after Quimpo 1983 JWRPM
                e_s = (y/b).^(-1/a);
            otherwise   % Now write Y as a function of e
                y_s = b * (e).^(-a);
        end
    case 'quimpo' % Quimpo - 2 parameters
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % FDC function
                e_s = -1/b * log(y/a);
            otherwise   % Now write Y as function of e
                y_s = a * exp(-b * e);
        end
    case 'viola' % Viola et al. - 2 parameters: 
        % http://www.hydrol-earth-syst-sci.net/15/323/2011/hess-15-323-2011.pdf
        % Relative duration, D is equivalent to e!!, and D_w equal to e_0
        % Hence relative duration zero streamflow days is 1 - e_0, see Page
        % 325 of link
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % FDC function
                e_s = 1./((y/a).^(1/b) + 1);
            otherwise   % Now write Y as function of e
                y_s = a * (1./e - 1).^b;
        end
    case 'genuchten-2' % van Genuchten - 2 parameters
        % Define c
        c = 1 - 1/b; 
        if any([a b c] < 0)
            [e_s,y_s] = deal(zeros(size(e)));
        else
            % Calibrate in e or y-space?
            switch strcmp(FDCPar.form,'e')
                case 1      % After van Genuchten retention function (1980)
                    e_s = (1 + (a * y).^b).^(-c);
                otherwise   % Now write Y as function of e
                    y_s = (1/a) * ((e.^(-1/c) - 1).^(1/b));
            end
        end
    case 'kosugi-2' % Kosugi - 2 parameters
        % Define c
        c = 0;
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % FDC function - after Kosugi (1994, 1996)
                e_s = ones(size(y));
                % Now only calculate if (c - y)/(c - a) > 0
                dummy = (y - c)./(a - c); id = find(dummy > 0);
                % Update (for dummy < 0 --> log gives imaginary numbers
                e_s(id) = 1/2 * erfc(log((y(id) - c) ./ (a - c)) / ...
                    (sqrt(2) * b));
            otherwise
                y_s = c + (a - c) * exp(erfcinv(2 * e) * sqrt(2) * b);
        end
    case 'exponential' % exponential - 2 parameters: http://goo.gl/0tfEJw
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % FDC function - after Quimpo 1983 JWRPM:
                e_s = (1/a * log(y/b));
            otherwise   % Now write Y as a function of e
                y_s = b * exp(a * e);
        end
    case 'gamma' % Gamma dist - 2 parameters: http://goo.gl/PzWELa
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % FDC function - After Leboutillier 1993 WRR
                e_s = 1 - cdf('Gamma',y,a,b);
            otherwise   % Now write Y as a function of e
                y_s = icdf('Gamma',1 - e(e <= e_0)/e_0,a,b);
                y_s(e > e_0) = 0;
        end
       
    %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                       3 PARAMETER EXPRESSIONS                       %
    %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    case 'lognormal-3' % lognormal - 3 parameters: http://goo.gl/PzWELa
        c = x(3);
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % FDC function - After Leboutillier 1993 WRR
                e_s = ones(size(y)); 
                % Now find values of Y > c,
                id = find(y > c);
                % Update e
                e_s(id) = 1 - 1/2 * erfc(-(log(y(id) - c) - a) / ...
                    (sqrt(2) * b));
            otherwise   % Now write Y as a function of e
                y_s = c + exp(a - sqrt(2) * b * erfcinv(2 * (1 - e)));
        end
    case 'pareto' % Generalized Pareto dist: 3 prmtrs: http://goo.gl/YRBASP
        c = x(3);
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % FDC function - after Castellarin 2004 WRR:
                theta = max((y - a)/b,0);
                % if c ~= 0
                cdf_p = 1 - (1 + c * theta).^(-1/c);
                % Now check
                if c < 0, id =  theta >- 1/c ; cdf_p(id) = 1; end
                % Calculate exceedance probability
                e_s = 1 - cdf_p;
            otherwise
                % 1 - e/e_0 = 1 - ( 1 + c * theta ) .^ ( -1 / c )
                % (e/e_0).^-c = ( 1 + c * theta )
                % 1/c * ( (e/e_0).^-c  - 1 ) = theta
                % ( y - a ) / b = 1/c * ( (e/e_0).^-c  - 1 )
                % y = a + b/c * ( (e/e_0).^-c  - 1 )

                % Now write Y as a function of e
                y_s = a + b/c * (e.^(-c) - 1);
                % Correct Y if e > e_0
                id = e > 1;
                % Update Y
                y_s(id) = 0;
        end
    case 'gev' % GEV distribution - 3 parameters
        % Define c
        c = x(3);
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % Check value of k
                if c > 0
                    % Need to check first ( k ( = x(3) ) > 0 )
                    y_check = a - b/c; id = find(y < y_check); value = 0;
                else
                    y_check = a + b/(-c); id = find(y > y_check); value = 1;
                    % Set id values equal to zero
                end
                % Calculate cdf
                cdf_gev = exp(- (1 + c * ((y - a)/b)).^(-1/c) );
                % Set cdf values
                cdf_gev(id) = value;
                % How calculate exceedance probability
                e_s = 1 - cdf_gev;
            otherwise
                % cdf = 1 - e/e_0; 
                % cdf = exp(- (1 + c * ((y - a)/b)).^(-1/c) )
                % 1 - e/e_0 = exp(- (1 + c * ((y - a)/b)).^(-1/c))
                % -log(1 - e/e_0) = (1 + c * ((y - a)/b)).^(-1/c)
                % (- log(1-e/e_0)).^(-c) = 1 + c * (y - a)/b
                % 1/c * ((- log(1-e/e_0) ).^(-c) - 1) = (y - a)/b
                % b/c * ((- log(1-e/e_0) ).^(-c) - 1) = (y - a)
                % y = a + b/c * ((- log(1-e/e_0)).^(-c) - 1)
            
                % Check value of k
                if c > 0
                    % Need to check first ( k ( = x(3) ) > 0 )
                    y_check = a - b/c;
                    % Now write Y as function of e
                    y_s = a + b/c * ((-log(1 - e)).^(-c) - 1);
                    % Now check
                    id = y < y_check ; y_s(id) = a - b/c;
                else
                    y_check = a + b/(-c);
                    % Now write Y as function of e
                    y_s = a + b/c * ((-log(1 - e)).^(-c) - 1);
                    % Now check
                    id = y > y_check ; y_s(id) = inf;
                end
        end
    case 'franchini' % Franchini and Suppo - 3 parameters
        c = x(3);
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % FDC function
                e_s = ones(size(y));
                % Now calculate dummy variable
                dummy = (y - b)/a; 
                % Find dummy > idx
                id = find(dummy > 0);
                % How calculate exceedance probability
                e_s(id) = 1 - (dummy(id)).^(1/c);
            otherwise
                % Now write Y as function of e
                y_s = b + a * (1 - e).^(c);
        end
    case 'genuchten-3' % van Genuchten - 3 parameters
        % Define c
        c = x(3);
        if any([a b c] < 0)
            [e_s,y_s] = deal(zeros(size(e)));
        else
            % Calibrate in e or y-space?
            switch strcmp(FDCPar.form,'e')
                case 1      % After van Genuchten retention function (1980)
                    e_s = (1 + (a * y).^b).^(-c);
                otherwise   % Now write Y as function of e
                    y_s = (1/a) * ((e.^(-1/c) - 1).^(1/b));
            end
        end
    case 'kosugi-3' % Kosugi - 3 parameters
        % Define c
        c = x(3);
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % FDC function - after Kosugi (1994, 1996)
                e_s = ones(size(y));
                % Now only calculate if (c - y )/(c - a) > 0
                dummy = (y - c)./(a - c); id = find(dummy > 0);
                % Update (for dummy < 0 --> log gives imaginary numbers
                e_s(id) = 1/2 * erfc(log(dummy(id))/(sqrt(2) * b));
            otherwise
                y_s = c + (a - c) * exp(erfcinv(2 * e) * sqrt(2) * b);
        end

    %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                   ADDITIONAL MODELS NOT USED HEREIN                 %
    %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        
    case 'weibull' % weibull dist - 2 parameters: 
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % FDC function - Weibull distribution
                %e = 1 - ( 1 - exp(- (y/a)^b)); -> Weibull
                e_s = exp(-(y/a).^b);
            otherwise
                % log(e) = - ( y / a)^{b}
                % -log(e)^(1/b) = ( y / a);
                y_s = -a * (log(e)).^(1/b);
        end
    case 'rayleigh' % rayleigh dist - 1 parameter 
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % FDC function - Rayleigh distribution
                %e = 1 - ( 1 - exp( - 1/2 * ( y.^2 / a^2 ) ) ); -> Weibull
                e_s = exp(- 1/2 * (y.^2/a^2));
            otherwise
                % log(e) = - 1/2 * ( y.^2 / a^2)
                % -2 * log(e) = ( y.^2 / a.^2 );
                y_s = sqrt(-2 * a^2 * log(e));
        end
    case 'non_central_F' % Noncentral F - 3 parameters
        c = x(3);
        % Calibrate in e or y-space?
        switch strcmp(FDCPar.form,'e')
            case 1      % FDC function - Rayleigh distribution
                %e = 1 - (1-exp(-1/2*(y.^2/a^2))); -> Weibull
                e_s = 1 - ncfcdf(y,a,b,c);
            otherwise
                error('Not available');
                % log(e) = - 1/2*(y.^2/a^2)
                % -2*log(e) = (y.^2/a.^2);
                % Y = sqrt(-2*a^2*log(e));
        end
    %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                 END ADDITIONAL MODELS NOT USED HEREIN               %
    %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
end

% Compute residuals
if strcmp(FDCPar.form,'e')
    % Need to multiply e with e_0 (e_0: exceedance probability zero flows)
    res = e - e_s;  % sim_data = real(e_s); 
else
    res = y - y_s;  % sim_data = real(y_s); 
end

switch returnarg
    case 0 % Exceedance probabilities
        switch strcmp(FDCPar.form,'e')
            case 1
                varargout(1) = {e_s};
            otherwise
                varargout(1) = {y_s};
        end
    case 1 % L2 objective function
        varargout(1) = {res'*res};
    case 2 % Residual vector
        varargout(1) = {res};
    case 3  % L1 objecitve function
        varargout(1) = {norm(res,1)}; 
    otherwise
        error('FDCFIT_models:Unknown model')
end