function [X,RMSE,it] = cmaes_code(FDCPar,e,y,options)
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
% (µ/µ_w,lambda)-CMA-ES algorithm                                         %
%                                                                         %
% SYNOPSIS: [X,RMSE,it] = cmaes_code(FDCPar,e,y,options)                  %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Dec. 2014                               %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

d = FDCPar.d;         % number of objective variables/problem dimension
x_mean = zeros(d,1);  % rand(N,1);   % objective variables initial point
sigma = 0.3;          % coordinate wise standard deviation (step size)
stopfitness = 1e-10;  % stop if fitness < stopfitness (minimization)
stopeval = 1e4*d^2;   % stop after stopeval number of function evaluations

% Strategy parameter setting: Selection
lambda = options.P;
mu = lambda/2;                  % # parents/points for recombination
w = log(mu+1/2)-log(1:mu)';     % weighted recombination
mu = floor(mu);
w = w/sum(w);                   % normalize recombination weights array
mueff = sum(w)^2/sum(w.^2);     % variance-effectiveness of sum w_i x_i

% Strategy parameter setting: Adaptation
cc = (4+mueff/d) / ...          % time constant for cumulation for C
    (d+4 + 2*mueff/d);  
cs = (mueff+2) / (d+mueff+5);   % t-const for cumulation for sigma control
c1 = 2 / ((d+1.3)^2+mueff);     % learning rate for rank-one update of C
cmu = min(1-c1, 2 * ...         % and for rank-mu update
    (mueff-2+1/mueff) / ...
    ((d+2)^2+mueff));  
damps = 1 + 2*max(0, ...        % damping for sigma
    sqrt((mueff-1)/ ...
    (d+1))-1) + cs;             % usually close to 1

% Initialize dynamic (internal) strategy parameters and constants
pc = zeros(d,1); ps = zeros(d,1);   % evolution paths for C and sigma
B = eye(d,d);                       % B defines the coordinate system
D = ones(d,1);                      % diagonal D defines the scaling
C = B * diag(D.^2) * B';            % covariance matrix C
invsqrtC = B * diag(D.^-1) * B';    % C^-1/2
eigeneval = 0;                      % track update of B and D
chiN = d^0.5*(1-1/(4*d)+1/(21*d^2));% expectation of
                                    %   ||N(0,I)|| == norm(randn(N,1))
child = nan(FDCPar.d,lambda); OFc = nan(1,lambda);
it = 0; count = 0; % twenty lines of interesting code below

while it < stopeval && (it < options.MaxFunEvals)
    % Generate and evaluate lambda offspring
    for k = 1:lambda
        % m + sig * Normal(0,C)
        child(:,k) = x_mean + sigma * B * (D .* randn(d,1));
        OFc(k) = FDCFIT_functions(child(:,k),FDCPar,e,y,1);
        % Update iteration counter
        it = it + 1;
    end
    % Sort by fitness and compute weighted mean into xmean
    [OFc, id_s] = sort(OFc); % minimization
    x_old = x_mean;
    x_mean = child(:,id_s(1:mu))*w;   % recombination, new mean value
    % Cumulation: Update evolution paths
    ps = (1-cs)*ps ...
        + sqrt(cs*(2-cs)*mueff) * ... 
        invsqrtC * (x_mean-x_old) / sigma;
    hsig = norm(ps)/sqrt(1-(1-cs)^ ...
        (2*it/lambda))/chiN < 1.4 + 2/(d+1);
    pc = (1-cc)*pc + hsig * ...
        sqrt(cc*(2-cc)*mueff) * ...
        (x_mean-x_old) / sigma;
    tmp = (1/sigma) * (child(:,id_s(1:mu)) ...
        - repmat(x_old,1,mu));
    C = (1-c1-cmu) * C ...              % Adapt covariance matrix C
        + c1 * (pc*pc' ...              % plus rank one update
        + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
        + cmu * tmp * diag(w) * tmp';   % plus rank mu update
    sigma = sigma * exp((cs/damps)* ... % Adapt step size sigma
        (norm(ps)/chiN - 1));
    
    % Decomposition of C into B*diag(D.^2)*B' (diagonalization)
    if it - eigeneval > ...             % to achieve O(N^2)
            lambda/(c1+cmu)/d/10  
        eigeneval = it;
        C = triu(C) + triu(C,1)';       % enforce symmetry
        [B,D] = eig(C);                 % eigen decomposition, 
                                        % B  are normalized eigenvectors
        D = sqrt(diag(D));              % D is vector standard devs. now
        invsqrtC = B * diag(D.^-1) * B';
    end

    if (it > 1)                         % Print progress
        fprintf(1,repmat('\b',1,count)); 
        count = fprintf(['FDCFIT CALCULATING: ' ...
            '%3.2f %% done of maximum of %u ' ...
            'trials (= options.MaxFunEvals) with %s'],...
            100*(it/options.MaxFunEvals), ...
            options.MaxFunEvals,'CMA-ES algorithm');
    end
    
    % Break, if fitness is good enough or condition exceeds 1e14, 
    % better termination methods are advisable
    if OFc(1) <= stopfitness ...
            || max(D) > 1e7 * min(D) ...
            || (max(max(child,[],2) - ...
            min(child,[],2))) < options.TolX ...
            || (OFc(options.P) - OFc(1)) < options.TolFun
            fprintf('\n'); fprintf(['FDCFIT CALCULATING: CMAES has ' ...
                'converged after %u generations\n'],it);
        break
    end
    
end % while, end generation loop

% Return parameter values
X = child'; RMSE = sqrt(OFc'/FDCPar.n);
% Repeat iter (trick for later)
it = it * ones(options.P,1);

end

% ---------------------------------------------------------------
% % function f=frosenbrock(x)
% %     if size(x,1) < 2 error('dimension must be greater one'); end
% %     f = 100*sum((x(1:end-1).^2 - x(2:end)).^2) + sum((x(1:end-1)-1).^2);