function [temp,P,rho] = atmosphere(H,Temp0,P0,tempr0,Pr0,H0,kTemp,R)
    if H < 36152 && H >= 0
        temp = 59 - 0.00356*H;
        tempr = temp;
        temp = tempr+(Temp0-tempr0)*exp((H0-H)/kTemp);
        P = 2116*((temp+459.67)/518.6)^5.256;
        Pr = P;
        P = Pr+(P0-Pr0)*exp((H0-H)/kTemp);
    elseif H >= 36152 && H < 82854
        temp = -70;
        tempr = temp;
        temp = tempr+(Temp0-tempr0)*exp((H0-H)/kTemp);
        P = 473.1*exp(1)^(1.73-0.000048*H);
        Pr = P;
        P = Pr+(P0-Pr0)*exp((H0-H)/kTemp);
    else
        temp = -205.05+0.00164*H;
        tempr = temp;
        temp = tempr+(Temp0-tempr0)*exp((H0-H)/kTemp);
        P = 51.97*((temp+459.67)/389.98)^-11.388;
        Pr = P;
        P = Pr+(P0-Pr0)*exp((H0-H)/kTemp);
    end
    rho = P/(R*(temp+459.67));
end