function J = FDCFIT_jac(FDCPar,x,y)
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
% Calculates the Jacobian of the VG or Kosugi FDC formulations using      %
% analytical functions of the partial parameter sensitivities (zero cost) %
%                                                                         %
% SYNOPSIS: J = FDCFIT_jac(FDCPar,x,Y)                                    %
%  where                                                                  %
%   FDCPar    [input] Structure with settings for FDCFIT                  %
%    .d           # unknown parameters of parametric FDC expression       %
%    .form        Forward or inverse formulation of the FDC               %
%      = 'e'      Forward expression: e_s = f(y)                          %
%      = 'y'      Backward/inverse expression: y_s = f(e)                 %
%   x         [input] 1xd vector parameter values parametric expression   %
%   y         [input] nx1 vector of empirical (= measured) discharges     %
%   J         [outpt] nxd Jacobian matrix of the FDC expression           %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Dec. 2014                               %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Now determine which model to use
switch FDCPar.model_class
    case {'vg'}     % van Genuchten
        % And which of the three formulations
        switch FDCPar.model_formulation
            case '2'    % two-parameter van Genuchten soil water 
                        % retention function (1980)
                % Define a and b
                a = x(1); b = x(2); c = 1 - 1/b;
                % p_pred = ( 1 + ( a * Y  ) .^b ) .^ ( -c );
                
                % dGdalpha = diff('( 1 + ( a * Y  ) ^b ) ^ ( -c )','a')
                J(:,1) = (c*b*y.*(a*y).^(b - 1))./((a*y).^b + 1).^(c + 1);
                % dGdn = diff('( 1 + ( a * Y  ) ^b ) ^ ( -c )','b')
                J(:,2) = (c*log(a*y).*(a*y).^b)./((a*y).^b + 1).^(c + 1);
            case '3'    % three-parameter van Genuchten soil water 
                        % retention function (1980)
                % Define a and b
                a = x(1); b = x(2); c = x(3);
                % p_pred = ( 1 + ( a * Y  ) .^b ) .^ ( -c );
                
                % dGdalpha = diff('( 1 + ( a * Y  ) ^b ) ^ ( -c )','a')
                J(:,1) = (c*b*y.*(a*y).^(b - 1))./((a*y).^b + 1).^(c + 1);
                % dGdn = diff('( 1 + ( a * Y  ) ^b ) ^ ( -c )','b')
                J(:,2) = (c*log(a*y).*(a*y).^b)./((a*y).^b + 1).^(c + 1);
                % dGdm = diff('( 1 + ( a * Y  ) ^b ) ^ ( -c )','c')
                J(:,3) = log((a*y).^b + 1)./((a*y).^b + 1).^c;
            case '5'    % five-parameter mixture formulation of two VG 
                        % functions for bimodal pore size distributions 
                        % (Durner, 1990)
                % Define a and b
                a = x(1); b = x(2); m_1 = 1 - 1/b; c = x(3); d = x(4); 
                m_2 = 1 - 1/d; w_1 = x(5); w_2 = 1 - w_1;
                % p_pred = w_1 * ( 1 + ( a * Y  ) .^b ) .^ ( -m_1 ) + ... 
                %     w_2 * ( 1 + ( c * Y  ) .^d ) .^ ( -m_2 );
                
                % dGdalpha_1 = diff(['w_1 * ( 1 + ( a * Y  )^b )^ ...' ...
                %     '( -m_1 ) + w_2 * ( 1 + ( c * Y  ) ^d ) ^ ' ...
                %     '( -m_2 )'],'a')
                J(:,1) = (m_1*b*y*w_1.*(a*y).^(b - 1))./((a*y).^b ...
                    + 1).^(m_1 + 1);
                % dGdn_1 = diff(['w_1 * ( 1 + ( a * Y  ) ^b ) ^ ' ...
                %     '( -m_1 ) + w_2 * ( 1 + ( c * Y  ) ^d ) ^ ' ...
                %     '( -m_2 )'],'b')
                J(:,2) = (m_1*w_1*log(a*y).*(a*y).^b)./((a*y).^b ...
                    + 1).^(m_1 + 1);
                % dGdalpha_2 = diff(['w_1 * ( 1 + ( a * Y  ) ^b ) ^' ...
                %     ' ( -m_1 ) + w_2 * ( 1 + ( c * Y  ) ^d ) ^' ...
                %     ' ( -m_2 )'],'c')
                J(:,3) = (m_2*d*y*w_2.*(c*y).^(d - 1))./((c*y).^d ...
                    + 1).^(m_2 + 1);
                % dGdn_2 = diff(['w_1 * ( 1 + ( a * Y  ) ^b ) ^' ...
                %     '( -m_1 ) + w_2 * ( 1 + ( c * Y  ) ^d ) ^' ...
                %     '( -m_2 )'],'d')
                J(:,4) = (m_2*w_2*log(c*y).*(c*y).^d)./((c*y).^d ...
                    + 1).^(m_2 + 1);
                % dGdw_1 = diff(['w_1 * ( 1 + ( a * Y  ) ^b ) ^ ' ...
                %     '( -m_1 ) + w_2 * ( 1 + ( c * Y  ) ^d ) ^ ' ...
                %     '( -m_2 )'],'w_1')
                J(:,5) = 1./((a*y).^b + 1).^m_1;
            otherwise
                disp('do not know this paticular formulation of VG');
        end
    case {'k'}      % Kosugi
        switch FDCPar.model_formulation
            case '2'    % two-parameter Kosugi lognormal pore size 
                        % distribution model (1994, 1996)
                a = x(1); b = x(2);
                % p_pred = 1/2 * erfc ( log ( Y / a ) / ( sqrt(2) * b ) );
                
                % dGda = diff(['1/2*erfc ( log ( Y / a ) / ' ...
                %     '( sqrt(2) * b ) )'],'a');
                J(:,1) = - 2^(1/2) ./ (2*pi^(1/2)*a*b* ...
                    exp(log(y/a).^2/(2*b^2)));
                % dGdb = diff(['1/2*erfc ( log ( Y / a )/(sqrt(2) * b ' ...
                %     ') )'],'b');
                J(:,2) = -(2^(1/2)*log(y/a))./(2*pi^(1/2)* ...
                    b^2*exp(log(y/a).^2/(2*b^2)));
            case '3'    % three-parameter Kosugi lognormal pore size 
                        % distribution model (1994, 1996)
                a = x(1); b = x(2); c = x(3);
                % p_pred = 1/2 * erfc ( log ( ( Y - c ) / ( a - c ) ) ...
                %     / ( sqrt(2) * b ) );
                J = zeros ( numel ( y ) , 3 );
                % Now determine
                dummy = ( y  - c ) ./ ( a - c ); id = find ( dummy > 0 );
                
                % dGda = diff(['1/2*erfc ( log ( ( Y - c ) / ( a - c ' ...
                %     ') ) / ( sqrt(2) * b ) )'],'a');
                J(id,1) = -2.^(1/2)./(2*pi.^(1/2)*b*exp(log(-(c - ...
                    y(id))./(a - c)).^2./(2*b^2))*(a - c));
                % dGdb = diff(['1/2*erfc ( log ( ( Y - c ) / ( a - c' ...
                %     ' ) ) / ( sqrt(2) * b ) )'],'b');
                J(id,2) = -(2.^(1/2)*log(-(c - y(id))./(a - ...
                    c)))./(2*pi.^(1/2)*b^2*exp(log(-(c - y(id))./ ...
                    (a - c)).^2./(2*b^2)));
                % dGdc = diff(['1/2*erfc ( log ( ( Y - c ) /(a - c)) ' ...
                %     '/ ( sqrt(2) * b ) )'],'c');
                J(id,3) = (2.^(1/2)*((c - y(id))./(a - c).^2 + ...
                    1./(a - c))*(a - c))./(2*pi.^(1/2)*b*exp( ...
                    log(-(c - y(id))./(a - c)).^2./(2*b^2)).*(c - y(id)));
            case '5'    % five-parameter mixture formulation two Kosugi WRF 
                        % functions for bimodal pore size distributions 
                        % (Durner, 1990)
                a = x(1); b = x(2); psi_c1 = 0; c = x(3); d = x(4); 
                w_1 = x(5); w_2 = 1 - w_1; 
                % psi_c2 = 0 ( replace psi_c1 of 2nd compnnt with psi_c2 )
                % p_pred = w_1 * 1/2 * erfc ( log ( ( Y - psi_c1 ) / ...
                %     ( a - psi_c1 ) ) / ( sqrt(2) * b ) ) + ...
                %     w_2 * 1/2 * erfc ( log ( ( Y - psi_c1 ) / ...
                %     ( c - psi_c1 ) ) / ( sqrt(2) * d ) );
                % dGda = diff(['w_1 * 1/2 * erfc ( log ((Y - psi_c1) ' ...
                %     '/ ( a - psi_c1 ) ) / ( sqrt(2) * b ) ) + ' ...
                %     'w_2 * 1/2 * erfc ( log ( ( Y - psi_c1 ) / ' ...
                %     '( c - psi_c1 ) ) / ( sqrt(2) * d ) )'] ,'a');
                J(:,1) = -(2^(1/2)*w_1)./(2*pi^(1/2)*b*exp(log(- ...
                    (psi_c1 - y)./(a - psi_c1)).^2./(2*b^2))*(a - psi_c1));
                % dGdb = diff(['w_1 * 1/2 * erfc ( log ((Y - psi_c1) ' ...
                %     '/ ( a - psi_c1 ) ) / ( sqrt(2) * b ) ) + ' ...
                %     'w_2 * 1/2 * erfc ( log ( ( Y - psi_c1 ) / ' ...
                %     '( c - psi_c1 ) ) / ( sqrt(2) * d ) )'] ,'b');
                J(:,2) = -(2^(1/2)*w_1*log(-(psi_c1 - y)./(a - ...
                    psi_c1)))./(2*pi^(1/2)*b^2*exp(log(-(psi_c1 - y)./ ...
                    (a - psi_c1)).^2./(2*b^2)));
                % dGdc = diff(['w_1 * 1/2 * erfc (log((Y - psi_c1) ' ...
                %     '/ ( a - psi_c1 ) ) / ( sqrt(2) * b ) ) + ' ...
                %     'w_2 * 1/2 * erfc ( log ( ( Y - psi_c1 ) / ' ...
                %     '( c - psi_c1 ) ) / ( sqrt(2) * d ) )'] ,'c');
                J(:,3) = -(2^(1/2)*w_2)./(2*pi^(1/2)*d*exp( ...
                    log(-(psi_c1 - y)./(c - psi_c1)).^2./(2*d^2))* ...
                    (c - psi_c1));
                % dGdd = diff(['w_1 * 1/2 * erfc ( log ( (Y - psi_c1) ' ...
                %     '/ ( a - psi_c1 ) ) / ( sqrt(2) * b ) ) + ' ...
                %     'w_2 * 1/2 * erfc ( log ( ( Y - psi_c1 ) / ' ...
                %     '( c - psi_c1 ) ) / ( sqrt(2) * d ) )'] ,'d');
                J(:,4) = -(2^(1/2)*w_2*log(-(psi_c1 - y)./(c - ...
                    psi_c1)))./(2*pi^(1/2)*d^2*exp(log(-(psi_c1 - ...
                    y)./(c - psi_c1)).^2./(2*d^2)));
                % dGdw_1 = diff(['w_1 * 1/2 * erfc ( log((Y - psi_c1 ' ...
                %  ') / ( a - psi_c1 ) ) / ( sqrt(2) * b ) ) + w_2 * ' ...
                %  '1/2 * erfc ( log ( ( Y - psi_c1 ) / ( c - psi_c1 ' ...
                %  ') ) / ( sqrt(2) * d ) )'] ,'w_1');
                J(:,5) = -erfc((2^(1/2)*log(-(psi_c1 - y)./(a - ...
                    psi_c1)))./(2*b))/2;
            otherwise
                disp('do not know this particular formulation of Kosugi');
        end
    otherwise
        disp('do not know this water retention function');
end

end
