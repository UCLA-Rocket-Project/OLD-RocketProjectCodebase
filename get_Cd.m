function [Cd,Mach] = get_Cd(time,burn_time,velocity,Temp,RAS)
%Finds Cd and Mach number using RASAero drag file, interpolation.

%Find Speed of sound,
gamma = 1.401;
R = 1716.49; %Air gas constant (ft lbf/slug R)
Temp = Temp + 491.67; %Convert to Rankine
v_Sound = sqrt(gamma*R*Temp); %in ft/s
%Find Mach Number

Mach = velocity/v_Sound;

M_l = ceil(100*Mach)/100;
M_h = M_l + 0.01;
if Mach == 0 %Causes indexing error -> use backwards interpolation
    M_h = M_h + 0.01;
    M_l = M_l + 0.01;
end

if time > burn_time
    CD_h = RAS(round(M_h*100), 4);
    CD_l = RAS(round(M_l*100), 4);
else
    CD_h = RAS(round(M_h*100), 5);
    CD_l = RAS(round(M_l*100), 5);
end

m = (CD_h-CD_l)/(M_h-M_l); %interpolation slope
b = CD_h - m*M_h; %interpolation intercept

Cd = m*Mach + b;

end
