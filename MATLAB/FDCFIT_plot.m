function FDCFIT_plot(FDCPar,e,y,method,map,RMSE_map,options,str, ...
    units,math_str)
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
% Plots the results of the flow duration curve formulation                %
%                                                                         %
% SYNOPSIS: FDCFIT_plot(FDCPar,e,y,method,map,RMSE_map,options,str, ...   %
%               units,math_str)                                           %
%  where                                                                  %
%   FDCPar    [input] Structure with settings for FDCFIT                  %
%   e         [input] nx1 vector of empirical exceedance probabilities    %
%   y         [input] nx1 vector of empirical (= measured) discharges     %
%   method    [input] Name (string) of parameter estimation method        %
%   map       [input] 1xd vector optmzed parameters parametric FDC expr.  %
%   RMSE_map  [input] scalar with minimized objective function            %
%   options   [input] Optional structure for screen output                %
%   str       [input] string with model equation FDC expression in latex  %
%   units     [input] string with model equation FDC expression in latex  %
%   math_str  [input] string with model equation FDC expression in latex  %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Dec. 2014                               %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

close all hidden;

set(0,'defaultTextInterpreter','latex'); % trying to set the default

% Define name of program
n_program = 'FDCFIT';
% Define name of figures file
file_name = [n_program,'_figures.pdf'];
% Determine screen size
scr_z = get(0,'ScreenSize');
% Multiplier, x and y: axis
x_mult = scr_z(3)/1920; y_mult = scr_z(4)/1080;
% Multiplier, text
t_mult = min(x_mult,y_mult);
% Define fontsize for figures
fontsize_titlepage = 30 * t_mult;
fontsize_labels = 21 * t_mult;
fontsize_axis = 17 * t_mult;
fontsize_legend = 15 * t_mult;
fontsize_marker = 20 * t_mult;
linewidth = 2 * t_mult;

% -------------------------------------------------------------------------

% Delete the warning file
dos('del FDCFIT_output.txt');

% open an output file with results of formulation
fid = fopen('FDCFIT_output.txt','w');
fprintf(fid,['------------------- FDCFIT output file -----------------' ...
    '---\n']);
fprintf(fid,'\n');
if strcmp(FDCPar.form,'e')
    fprintf(fid,['→ Results for model %s (Check Table 2 of ' ...
        'manual) \n'],upper(char(FDCPar.model)));
else
    fprintf(fid,['→ Results for model %s (Check Table 1 of ' ...
        'manual) \n'],upper(char(FDCPar.model)));
end
%fprintf(fid,'\n');
switch FDCPar.form
    case 'y'
        fprintf(fid,['→ Flow duration curve minimizes' ...
            ' discharge residuals\n']);
    case 'e'
        fprintf(fid,['→ Flow duration curve minimizes' ...
            ' exceedance probability residuals\n']);

end
%fprintf(fid,'\n');
fprintf(fid,'→ Parameter estimation method: %s algorithm\n', ...
    upper(char(method)));
fprintf(fid,'\n');
fprintf(fid,'==================================\n');
if numel(char(FDCPar.acronym)) == 1
    fprintf(fid,'Parameter      Value      Units \n');
elseif numel(char(FDCPar.acronym)) == 2
    fprintf(fid,'Parameter       Value      Units \n');
elseif numel(char(FDCPar.acronym)) == 3
    fprintf(fid,'Parameter        Value      Units \n');
end
fprintf(fid,'----------------------------------\n');
% Now print parameter values
for i = 1 : FDCPar.d
    % Create parameter name
    evalstr = strcat(str(i),FDCPar.acronym);
    % Now print to file
    if map(i) < 0
        fprintf(fid,' %s         %7.5f     %s\n',char(evalstr), ...
            map(i),char(units(i)));
    elseif ( map(i) >= 10 ) && ( map(i) < 100 )
        fprintf(fid,' %s          %7.4f     %s\n',char(evalstr), ...
            map(i),char(units(i)));
    elseif ( map(i) >= 100 ) && ( map(i) < 1000 )
        fprintf(fid,' %s          %7.3f     %s\n',char(evalstr), ...
            map(i),char(units(i)));
    elseif ( map(i) >= 1000 )
        fprintf(fid,' %s          %7.2f     %s\n',char(evalstr), ...
            map(i),char(units(i)));
    else
        fprintf(fid,' %s          %7.5f     %s\n',char(evalstr), ...
            map(i),char(units(i)));
    end
end
fprintf(fid,'----------------------------------\n');
if numel(char(FDCPar.acronym)) == 1
    if strcmp(FDCPar.form,'e')
        fprintf(fid,[' RMSE         %7.5f     (-)     ' ...
            '--> in exceedance probability space\n'],RMSE_map);
    else
        fprintf(fid,[' RMSE         %7.5f     (L/T)   ' ...
            '--> in discharge space\n'],RMSE_map);
    end
elseif numel(char(FDCPar.acronym)) == 2
    if strcmp(FDCPar.form,'e')
        fprintf(fid,[' RMSE          %7.5f     (-)     ' ...
            '--> in exceedance probability space\n'],RMSE_map);
    else
        fprintf(fid,[' RMSE          %7.5f     (L/T)   ' ...
            '--> in discharge space\n'],RMSE_map);
    end
elseif numel(char(FDCPar.acronym)) == 3
    if strcmp(FDCPar.form,'e')
        fprintf(fid,[' RMSE           %7.5f     (-)     ' ...
            '--> in exceedance probability space\n'],RMSE_map);
    else
        fprintf(fid,[' RMSE           %7.5f     (L/T)   ' ...
            '--> in discharge space\n'],RMSE_map);
    end
end
fprintf(fid,'==================================\n');
fprintf(fid,'\n');
% Write final line of warning file
fprintf(fid,['---------------- End of FDCFIT output file -------------' ...
    '---']);
% Close the warning file
fclose('all');
% Open the warning file
if (ispc || ismac), edit FDCFIT_output.txt, end

% ----------------------------------------------------------------------- %
%                 Now plot empty figure for PDF file                      %
% ----------------------------------------------------------------------- %

figure('units','normalized','outerposition',[0 0 1 1],'name','first page');
plot([],[],'ro'); axis([0 1 0 1]); set(gcf,'color','w');
set(gca,'XColor','w','YColor','w');
%title('Visual results of FDCFIT toolbox','fontsize',20,'interpreter','latex');
text(0.3*x_mult,0.6*y_mult,'Visual results of FDCFIT toolbox','fontsize', ...
    fontsize_titlepage,'interpreter','latex');
% Now about Tables
% text(0.3*x_mult,0.5*y_mult,'GUI Tables may not print well ', ...
%    'fontsize', fontsize_titlepage,'interpreter','latex');

% ----------------------------------------------------------------------- %

% Determine scaling
Y_max = floor(1.02 * max(y)); Y_min = max(1e-4,min(y));
X_min = -0.04; X_max = 1.04;
% Now maximize figure to screen
% set(gcf,'units','normalized','outerposition',[0 0 20 10]);
% figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
% determine screen size
% scrsz = get(groot,'ScreenSize');
% make figure
% figure('Position',[50 scrsz(4)/10 scrsz(3)-100 scrsz(4)/1.5]);
% axpos_old = [50 scrsz(4)/10 scrsz(3)-100 scrsz(4)/1.5];
fig_pos = [0.05 0.1 0.9 2/3]; figure('units','normalized', ...
    'position',fig_pos)
% Define axis
ax1 = axes('units','inches'); axpos1 = [10 1 7.5 6];
set(ax1,'position',axpos1);

% Now plot the model fit corresponding to the best formulation coefficients
if strcmp(FDCPar.form,'e')
    e_s = FDCFIT_functions(map,FDCPar,e,y,0);
else
    y_s = FDCFIT_functions(map,FDCPar,e,y,0);
end

%%%%%%%%%%%%%%%%%%%%%%% Logarithmic y-scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now plot the results --> first the measured FDC
h1 = semilogy(e,y,'r.'); set(h1,'markersize',fontsize_marker); hold on;
% Now plot the model fit corresponding to the best formulation coefficients
if strcmp(FDCPar.form,'e')
    semilogy(e_s,y,'b','linewidth',linewidth);
else
    semilogy(e,y_s,'b','linewidth',linewidth);
end
set(gca,'xtick',0.0:0.1:1.0,'xticklabel',{'0.0','0.1','0.2', ...
    '0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});
axis([X_min X_max Y_min Y_max]);
% Now add legend, axis labels, etc.
% [~,objh] = legend(evalstr_legend,'interpreter','latex');
%  legend boxoff; set(objh,'linewidth',3);
% try set(objh,'interpreter','latex','linewidth',3, ...
%     'markersize',30); catch, end;
% legend({'data',math_str},'interpreter','latex', ...
%     'location','northeast'); legend boxoff;

% ylabel('${\rm STREAMFLOW,} \hspace{1mm} \widetilde{\textbf{Y}} ...
%     \hspace{2mm} {\rm (L/T)}$','fontsize',16,'interpreter','latex');

% Increase fontsize to 16
set(gca,'fontsize',fontsize_axis,'tickdir','out');

% And scale curve
%if Y_min > 0
%axis([ X_min X_max 1e-3 Y_max]);
%else
%    axis([ X_min X_max 1e-3 Y_max]);
%end
% if Y_min > 0
%     % Set axes
%     axis([ 0 1 0 exp(log(Y_min) + 11/10 * (log(Y_max) - log(Y_min)))]);
% else
%     axis([ 0 1 0 exp(log(1e-2) + 11/10 * (log(Y_max) - log(1e-2)))]);
% end

% Add RMSE of the error of the fit to the figure
a = num2str(RMSE_map); idx_zeros = find ( a == num2str(0) ...
    | a == '.' ); ii = find ( diff ( idx_zeros ) > 1 );
if isempty(ii)
    last_zero = idx_zeros(end);
else
    last_zero = idx_zeros(ii);
end
if isempty(last_zero)
    % Just plot the first four values of a
    idx_plot = (1:min(find(a == '.') + 3,numel(a)));
else
    % Check whether we start with zero or not?
    idx_plot = (1:min( numel(a),last_zero + 3));
end

% Add legend
yr = ylim; y_pos1 = 10^(9.5/10 * (log10(yr(2)) - ...
    log10(yr(1))) + log10(yr(1)));
plot(0.07,y_pos1,'ro','MarkerFaceColor','red','markersize', ...
    fontsize_marker/2);
text(0.12,y_pos1,'Measured data','color','red','fontsize', ...
    fontsize_legend);
y_pos2 = 10^(8.5/10 * (log10(yr(2)) - ...
    log10(yr(1))) + log10(yr(1)));
line([0.04 0.10],[y_pos2 y_pos2],'color','blue','linewidth',2*linewidth)
text(0.12,y_pos2,math_str,"FontSize",fontsize_legend,"Color",'blue');
% print RMSE
y_pos3 = 10^(7.5/10 * (log10(yr(2)) - ...
    log10(yr(1))) + log10(yr(1)));
if strcmp(FDCPar.form,'e')
    evalstr = strcat('RMSE = ',{' '},a(idx_plot(1):idx_plot(end)), ...
        {' '},'(-)'); text(0.12,y_pos3, evalstr, ...
        'fontsize',fontsize_legend,'interpreter','latex');
else
    evalstr = strcat('RMSE = ',{' '},a(idx_plot(1):idx_plot(end)), ...
        {' '},'(L/T)'); text(0.12,...
        y_pos3,evalstr, ...
        'fontsize',fontsize_legend,'interpreter','latex');
end
xlh = xlabel(['${\rm Normalized \hspace{2mm} exceedance \hspace{2mm} ' ...
    'probability,} \hspace{2mm} e_{\rm n} (-)$'], ...
    'fontsize',fontsize_labels,'interpreter','latex');
% Change x-label position
xlh.Position(2) = xlh.Position(2) - 0.02*(log10(yr(2)) - log10(yr(1))) * yr(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Linear y-scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define axis
ax1 = axes('units','inches'); axpos1 = [1.5 1 7.5 6];
set(ax1,'position',axpos1);

% Now plot the results --> first the data
h2 = plot(e,y,'r.'); set(h2,'markersize',fontsize_marker); hold on;

% Now plot the model fit corresponding to the best formulation coefficients
if strcmp(FDCPar.form,'e')
    plot(e_s,y,'b','linewidth',linewidth);
else
    plot(e,y_s,'b','linewidth',linewidth);
end
set(gca,'xtick',0.0:0.1:1.0,'xticklabel',{'0.0','0.1','0.2', ...
    '0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});
axis([X_min X_max Y_min Y_max]);
% Increase fontsize to 16
set(gca,'fontsize',fontsize_axis,'tickdir','out');

% Add legend
yr = ylim; y_pos1 = 9.5/10 * (yr(2) - yr(1)) + yr(1);
plot(0.07,y_pos1,'ro','MarkerFaceColor','red','markersize', ...
    fontsize_marker/2);
text(0.12,y_pos1,'Measured data','color','red','fontsize', ...
    fontsize_legend);
y_pos2 = 8.5/10 * (yr(2) - yr(1)) + yr(1);
line([0.04 0.10],[y_pos2 y_pos2],'color','blue','linewidth',2*linewidth)
text(0.12,y_pos2,math_str,"FontSize",fontsize_legend,"Color",'blue');
% print RMSE
y_pos3 = 7.5/10 * (yr(2) - yr(1)) + yr(1);
if strcmp(FDCPar.form,'e')
    evalstr = strcat('RMSE = ',{' '},a(idx_plot(1):idx_plot(end)), ...
        {' '},'(-)'); text(0.12, ...
        y_pos3,evalstr, ...
        'fontsize',fontsize_legend,'interpreter','latex');
else
    evalstr = strcat('RMSE = ',{' '},a(idx_plot(1):idx_plot(end)), ...
        {' '},'(L/T)'); text(0.12, ...
        y_pos3,evalstr , ...
        'fontsize',fontsize_legend,'interpreter','latex');
end
xlh = xlabel(['${\rm Normalized \hspace{2mm} exceedance \hspace{2mm} ' ...
    'probability,} \hspace{2mm} e_{\rm n} (-)$'], ...
    'fontsize',fontsize_labels,'interpreter','latex');
xlh.Position(2) = xlh.Position(2) - 0.02*(yr(2) - yr(1));
ylh = ylabel(['${\rm Streamflow,} \hspace{1mm} y' ...
    '\hspace{2mm} {\rm (L/T)}$'], ...
    'fontsize',fontsize_labels,'interpreter','latex');
ylh.Position(1) = ylh.Position(1) - 0.02;


% Now add title
switch lower(options.type)
    case 'annual'
        FDC_type = lower(options.type);
    case 'day'
        FDC_type = 'daily';
    otherwise
        FDC_type = strcat(options.type,'ly');
end
switch FDCPar.form
    case 'e'
        % Fitted using exceeedance probabilities
        fit_str = strcat(['obtained by minimizing exceedance ' ...
            'probability residuals']);
    case 'y'
        % Fitted using exceeedance probabilities
        fit_str = strcat('obtained by minimizing discharge residuals');
end
evalstr = strcat(['FDCFIT Results V2.0: Comparison of empirical ' ...
    'and fitted'],{' '},FDC_type,{' '}, ...
    'flow duration curves',{' '},fit_str);
% Now add title
text(1.11,Y_min + 1.05*( Y_max - Y_min ),evalstr, ...
    'fontsize',18,'fontweight','bold','interpreter', ...
    'latex','HorizontalAlignment','center');
% Add test to right plot (print must appear at same height in figure)
uistack(h1,'down');

% NO PRINTING TO FILE
% % % Get figure handles
% % figHandles = flipud(findall(0,'Type','figure'));
% % for zz = 1:numel(figHandles)
% %     figure(figHandles(zz)); set(gcf,'color','w');
% %     switch zz
% %         case 1
% %             exportgraphics(figHandles(1),file_name);
% %         otherwise
% %             % Append to existing figure
% %             exportgraphics(figHandles(zz),file_name,'Append',true);
% %     end
% % end
% %
% % % Open PDF document
% % open(file_name);

% Print wait statement to the screen
fprintf('FDCFIT:done\n');

end
