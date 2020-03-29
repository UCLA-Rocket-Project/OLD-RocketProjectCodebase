function [noseCNa, bodyCNa, tailCNa] = getCNaNoseBodyTail(Ma, noseDiameter, noseLength, bodyDiameter, bodyLength, tailDiameter, tailLength, CNai0, aoai_i)
% Returns CNa value for any regime(transsonic, supersonic, even PG compressibility correction)

%Subsonic
    % Uses Prandtl-Glauert Compressibility Correction:
    % (http://cambridgerocket.sourceforge.net/AerodynamicCoefficients.pdf)
    if Ma < 0.7
        beta = sqrt(1-Ma^2); % PG correction
        noseCNa = CNai0(1)/beta;
        bodyCNa = CNai0(2)/beta;
        tailCNa = CNai0(3)/beta;
%Transonic
    % Currently no model:
    % Just plain Barrowman :(
    elseif Ma >= 0.7 && Ma < 1.2
        noseCNa = CNai0(1);
        bodyCNa = CNai0(2);
        tailCNa = CNai0(3);
       
%Supersonic
    elseif Ma > 1.2
        % Sim's exact method of characteristics (tabulated):
        % (https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19660028006.pdf)
            theta = rad2deg(atan((noseDiameter/2)/noseLength)); %approximate nose as cone, find cone angle
            m = 0.05525641*theta - 0.067948718; %two corresponding charts have theta = 5,10deg. So linearize for small changes in theta
            b = -0.044717949*theta + 2.03025641; 
            nosebodyCNa = m*Ma + b; %linear with Mach number until Ma ~ 4.0
            noseCNa = 0.2*Ma + 2.03; %approximate as convex parabolic ogive, linear with Mach number for 1.2<Ma<2.0
            bodyCNa = nosebodyCNa - noseCNa;
            tailCNa = CNai0(3); %Currently still no good model for boattail
        % Slender body approximation:
        % (Pasquale Manned Spacecraft Design Principles Section 7.6.2)
%             n = 1; CD_N = 0.53; aoai_avg = sum(aoai_i(1:3))/3;
%             noseSb = 0.25*pi*noseDiameter^2;
%             bodySb = 0.25*pi*bodyDiameter^2;
%             tailSb = 0.25*pi*tailDiameter^2;
%             S = bodySb;
%             noseSp = 0.5*noseDiameter*noseLength;
%             bodySp = bodyDiameter*bodyLength;
%             tailSp = 0.5*(bodyDiameter+tailDiameter)*tailLength;
% 
%             CN_eff = (tailSb/S)*sin(2*aoai_avg)*cos(aoai_avg/2)...
%                    + n*CD_N*((noseSp+bodySp+tailSp)/S)*(sin(aoai_avg))^2;
%             noseCN = (noseSb/S)*sin(2*aoai_i(1))*cos(aoai_i(1)/2)...
%                    + n*CD_N*(noseSp/S)*(sin(aoai_i(1)))^2;
%             nosebodyCN = (bodySb/S)*sin(sum(aoai_i(1:2)))*cos(sum(aoai_i(1:2))/4)... %%% Use better averaging method?
%                        + n*CD_N*((noseSp+bodySp)/S)*(sin(sum(aoai_i(1:2))/2))^2;
%             bodyCN = nosebodyCN - noseCN;
%             tailCN = CN_eff - nosebodyCN;
% 
%             noseCNa = noseCN/aoai_i(1);
%             bodyCNa = bodyCN/aoai_i(2);
%             tailCNa = tailCN/aoai_i(3);
    else
        error('Mach Number messed up: getCNaNoseBodyTail');
    end