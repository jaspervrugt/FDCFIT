function [X,RMSE,it] = de_code(X,FDCPar,e,y,Par_info,options)
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
% Differential Evolution Code                                             %
%                                                                         %
% SYNOPSIS: [X,RMSE,it] = de_code(X,FDCPar,e,y,Par_info,options)          %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Dec. 2014                               %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Replicate parameter boundaries
Par_info.min = repmat(Par_info.min,options.P,1); 
Par_info.max = repmat(Par_info.max,options.P,1);

% Start counter
it = options.P; tr = ceil(0.8*options.P);

% Calculate OF
OF = nan(options.P,1);
for i = 1:options.P
    OF(i,1) = FDCFIT_functions(X(i,1:FDCPar.d),FDCPar,e,y,1);
end

% Check for best solution
[~,id_s] = sort(OF); F = 0.6; G = 0.4; count = 0; converged = 0;

% Loop until converged
while ~converged && (it < options.MaxFunEvals)
    % Draw values to generate a, b and c from
    [~,draw] = sort(rand(options.P-1,options.P));
    % Create a, b and c
    r = draw(1:3,:)'; %F = 1/2 + 1/2*rand;
    % Create offspring assuming full crossover
    Z = X(1:options.P,1:FDCPar.d) + ...
        F * (X(r(1:options.P,2),1:FDCPar.d) - ...
        X(r(1:options.P,3),1:FDCPar.d)) + ...
        G * (X(id_s(1)*ones(options.P,1),1:FDCPar.d) - ...
        X(1:options.P,1:FDCPar.d) );
    % Now apply crossover
    id = rand(options.P,FDCPar.d) > options.CR;
    % Change those values of Z back to their parent value
    Z(id) = X(id);
    % Make sure that everything is in bound --> reflect or not
    id = find(Z < Par_info.min); Z(id) = 2 * Par_info.min(id) - Z(id);
    % Calculate OF for each child
    OFc = nan(options.P,1);
    for i = 1 : options.P
        OFc(i,1) = FDCFIT_functions(Z(i,1:FDCPar.d),FDCPar,e,y,1);
    end
    % Now accept or not
    id = OFc < OF;
    % Store new parents and objective function values
    X(id,1:FDCPar.d) = Z(id,1:FDCPar.d); OF(id) = OFc(id);
    % Check for convergence
    [OF_s,id_s] = sort(OF);
    % Now determine range of 95% of population (to avoid outlier problem)
    ii = id_s(1:tr);
    % Update counter
    it = it + options.P;
    % Print progress
    if (it > 1)
        fprintf(1, repmat('\b',1,count)); %delete line before
        count = fprintf(['FDCFIT CALCULATING: %3.2f %% done of ' ...
            'maximum of %u trials (= options.MaxFunEvals) with %s'],...
            100*(it/options.MaxFunEvals), ...
            options.MaxFunEvals,'Differential Evolution');
    end
    % Convergence achieved?
    if (max(abs(max(X(ii,:)) - min(X(ii,:)))) < options.TolX ) ...
            && ((OF_s(tr) - OF_s(1)) < options.TolFun)
        % Code has converged
        converged = 1; 
        fprintf('\n'); fprintf(['FDCFIT CALCULATING: Differential ' ...
            'Evolution has converged after %u generations\n'],it);
    end
    
end

% Calculate RMSE
RMSE = sqrt(OF/FDCPar.n);
% Repeat iter (trick for later)
it = it * ones(options.P,1);

end
