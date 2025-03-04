function [e_n,y_s,p_0] = calc_FDC(ID,type,col)
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
% Calculates the flow duration curve (FDC)                                %
%                                                                         %
% SYNOPSIS: [e_n,y_s,p_0] = calc_FDC(ID,type,col)                         %
%  where                                                                  %
%   ID        [input] String for watershed ID                             %
%   type      [input] Type of FDC derived from daily discharge data, y    %
%    = 'day'        Daily flow duration curve                             %
%    = 'week'       Weekly flow duration curve                            %
%    = 'month'      Monthly flow duration curve                           %
%    = 'annual'     Annual flow duration curve                            %
%   col       [input] column # of data with discharge data                %
%   e_n       [outpt] normalized exceedane probabilities for FDC type     %
%   y_s       [outpt] sorted discharge data for FDC type                  %
%   p_0       [outpt] probability of zero flows (of type)                 %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Dec. 2014                               %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if isnumeric(ID)
    error('calc_FDC: First input argument must be a matrix (.dly file)')
elseif ischar(ID)
    % Load the data
    data = load_data_dly(ID);
    % Take only zero or postive streamflow values ( should be all values)
    ii = data(:,col) >= 0; data = data(ii,:);
    % What is the size of daily data
    N = size(data,1);
end

% Now switch over possibilities of FDC
switch type
    case 'day'      % daily FDC
        % p_0 -> probability of zero flows
        ii = data(1:N,col) == 0; p_0 = sum(ii)/N;
        % Remove negative values
        Y = data(~ii,col);
    case 'week'     % weekly FDC
        % Set initial counter
        ct = 1; p = 0;
        % Average weekly data - take avg. of 7 days
        Y = nan(ceil(N/7),1);
        for i = 1:7:N
            % Extract data
            y = data(i:min(i + 6,N),col);
            % Remove negative values
            y = y(y > 0);
            if ~isempty(y)
                % Average streamflow values
                Y(ct,1) = mean(y);
                % Update counter
                ct = ct + 1;
            else
                p = p + 1;
            end
            % p_0 -> probability zero weekly flows
            p_0 = p/(ct-1);
        end

    case 'month'    % monthly FDC
        % Set initial counter
        ct = 1; p = 0;
        % Find months
        id_month = [1 ; find(diff(data(:,2)) ~= 0) + 1]; 
        % # months
        N_month = size(id_month,1);
        % Average weekly data
        Y = nan(N_month,1);
        for i = 2:N_month
            % Extract data
            y = data(id_month(i-1):id_month(i) - 1,col);
            % Remove negative values
            y = y(y > 0);
            if ~isempty(y)
                % Average streamflow values
                Y(ct,1) = mean(y);
                % Update counter
                ct = ct + 1;
            else
                p = p + 1;
            end
            % p_0 -> probability of zero monthly flows
            p_0 = p/(ct-1);
        end
    case 'annual'       % yearly FDC
        % Set initial counter
        ct = 1; p = 0;
        % Find months
        id_year = [1 ; find(diff(data(:,1)) ~= 0) + 1];
        % # years
        N_year = size(id_year,1);
        % Average yearly data
        Y = nan(N_year,1);
        for i = 2:size(id_year,1)
            % Extract data
            y = data(id_year(i-1):id_year(i) - 1,col);
            % Remove negative values
            y = y(y > 0);
            if ~isempty(y)
                % Average streamflow values
                Y(ct,1) = mean(y);         
                % Update counter
                ct = ct + 1;
            else
                p = p + 1;
            end
            % p_0 -> probability of zero yearly flows
            p_0 = p/(ct-1);
        end
    otherwise
        error(['FDCFIT:calc_FDC:wrong type ' ...
            '--> should be day/week/month/annual']);
end
% Remove nan due to initialization
Y = Y(~isnan(Y));
% How many values of Y?
N = size(Y,1);
% Now sort the discharge data in ascending order
y_s = sort(Y);
% Exceedance probabilities (based on N only)
e_n = 1 - ((1:N)' -0.5)./N;

end
