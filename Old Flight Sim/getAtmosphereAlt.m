function [rho,P,T] = getAtmosphereAlt(ry, h0, C1, C2, T1)
% Altered US Standard Atmosphere for range of launch atmospheric conditions
% Input//Output units: (slug/ft^3), (lbf/ft^2), (R//F)
    h1 = 36089; % (ft) Geopotential height station 1
    h2 = 65617; % (ft) Geopotential height station 2

    rE = 6378100/0.3048; % (ft) Radius of earth
    h = ry + h0; % (ft) Geometric elevation above sea level
    h = rE - ((rE^2)/(rE+h)); % (ft) Geopotential height (assumes spherical Earth for now)
    
    g = 32.174049; %(ft/s^2) gravity acceleration constant
    R = 1717; %(ft*lbf/slug*R) gas constant air
    
    if 0<=h && h<h1
        dTdh = -0.001981 * 9/5; % (R/ft)
        
        T = T1 - dTdh*(h1-h); % (R)
        P = exp(C1 - (g/(R*dTdh))*log(T)); % (lbf/ft^2) I think?
        rho = P/(R*T);
        
    elseif h1<=h && h<h2
        T = T1;
        P = exp(C2 - (g*h/(R*T)));
        rho = P/(R*T);
        
    else
        error('Exceeded max expected altitude');
    end
    T = T - 459.67; % (F)