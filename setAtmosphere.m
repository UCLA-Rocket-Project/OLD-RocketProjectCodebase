function [C1, C2, T1] = setAtmosphere(h0, T0, P0)
% Altered US Standard Atmosphere for range of launch atmospheric conditions
% Input//Output units: (slug/ft^3), (lbf/ft^2), (F//R)
    h1 = 36089; % (ft) Geopotential height station 1
    h2 = 65617; % (ft) Geopotential height station 2
    
    T0 = T0 + 459.67; % (R) Conversion to Rankine from Fahrenheit
    
    g = 32.174049; %(ft/s^2) gravity acceleration constant
    R = 1717; %(ft*lbf/slug*R) gas constant air
    
    %% First Integration Constant: Station 0 - 1
    % Assume h0<36089ft, which is pretty reasonable; FAR h<9000ft
    dTdh_1 = -0.001981 * 9/5; % (R/ft)
    C1 = log(P0) + (g/(R*dTdh_1))*log(T0); % Integration constant from analytical calcs
    T1 = T0 + dTdh_1*(h1-h0); % Temp calculation at station 1
    P1 = exp(C1 - (g/(R*dTdh_1))*log(T1)); % Pressure calculation at station 1
    
    %% Second Integration Constant: Station 1 - 2
    % In this zone T = T1 = constant, at boundary P = P1
    C2 = log(P1) + (g*h1/(R*T1)); % Use station 1 for boundary conditions
end