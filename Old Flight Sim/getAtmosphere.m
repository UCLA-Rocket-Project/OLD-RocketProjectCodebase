function [r,P,T] = getAtmosphere(h,h0)
%Atmospheric properties from NASA; uses altitude, pressure, temp
%Units: (slug/ft^3), (lbf/ft^2), (F)
%h0 is launch altitude
    alt = h+h0;

    if alt<36152
        T = 59-(0.00356*alt);
        P = 2116*(((T+459.7)/518.6)^5.256);
    end
    if alt>=36152 && alt<82345
        T = -70;
        P = 473.1*exp(1.73-(0.000048*alt));
    end
    if alt>=82345
        T = -205.05+(0.00164*alt);
        P = 51.97*(((T+459.7)/389.98)^-11.388);
    end
    r = P/(1718*(T+459.7));