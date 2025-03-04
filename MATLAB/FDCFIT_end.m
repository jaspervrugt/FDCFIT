function [map,RMSE_map,it_map,str] = FDCFIT_end(FDCPar,MAP,RMSE_MAP,it_MAP)
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
% Finalize return arguments, file writing, and setup calculation          %
%                                                                         %
% SYNOPSIS: [map,RMSE_map,it_map,str] = FDCFIT_end(FDCPar,MAP, ...        %
%               RMSE_MAP,it_MAP)                                          %
%  where                                                                  %
%   FDCPar    [input] Structure with settings for FDCFIT                  %
%   MAP       [input] Nxd matrix of N optimized parameter vectors         %
%   RMSE_MAP  [input] Nx1 vector of RMSE of N parameter vectors           %
%   it_MAP    [input] # iterations each of N multi-start trials           % 
%   map       [outpt] 1xd vector with overall best (map) parameter values %
%   RMSE_map  [outpt] scalar with RMSE of map parameter values            %
%   it_map    [outpt] # iterations of map                                 %
%   str       [outpt] string with parameter names                         %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Dec. 2014                               %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Find best values
[~,id] = min(RMSE_MAP);

% Extract lowest error and best parameters
map = MAP(id(1),1:FDCPar.d); RMSE_map = RMSE_MAP(id(1)); 
% # iterations of map parameter values
it_map = it_MAP(id(1));

% Define strig with names of parameters
str = {'a_','b_','c_'};

% open an output file with warnings
fid = fopen('warning_file.txt','a+','n');

% Write final line of warning file
fprintf(fid,'----------- End of MODELAVG warning file ----------\n');
% Close the warning file
fclose(fid);
% Now print to screen or not (not on unix/linux)
if ( ispc || ismac ), edit warning_file.txt, end

end
