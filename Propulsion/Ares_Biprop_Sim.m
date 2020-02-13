%% 
close all
% thermodynamic fluid properties and constants
mw_HE = 4;      %lb/lbmol
gam_HE = 1.67;
R = 10.73159;    % psi*ft3/lbmol-R
R_HE = R/mw_HE;  %psi*ft3/lb-R
rho_ox = 72.2;   %lb/ft3
rho_fuel = 50.6; % lb/ft3
g = 32.174;

% System Properties
m_ox = 28.16; %     lbm
m_fuel = 20.11; %   lbm

vol_ox = m_ox/rho_ox; % ft3
vol_fuel = m_fuel/rho_fuel; % ft3
vol = vol_ox+vol_fuel; %ft3
vol_t = 0.05787; % ft3
regset = 400; % regulator set pressure, psi
SPE = 0.013; % supply pressure effect

% geometries
L = 2; % ft
Ac = pi/4*4^2/144; % ft2
At = pi/4*1.5^2; % ft2
Vc = L*Ac; % ft3

cstar = 4132.8; % ft/s

% Flow Resistances
CdA_lox_inj = 0.0531; % in2
CdA_fuel_inj = 0.03703; % in2
CdA_lox_lines = 0.05175; % in2
CdA_fuel_lines = 0.05274; % in2
CdA_lox = 1/(sqrt(1/CdA_lox_inj^2+1/CdA_lox_lines^2)); %in2
CdA_fuel = 1/(sqrt(1/CdA_fuel_lines^2+1/CdA_fuel_inj^2));  %in2

% Sim Parameters
t_end = 1;
dt = 0.0001;
t = 0:dt:t_end;

% Stored Values
Pc_arr = zeros(1,length(t));
T_arr = zeros(1,length(t));
OF_arr = zeros(1,length(t));

% Initial Conditions
Pt = 4500; % psi
T_c = 300*1.8; % R
Pc = 50; % psi
OF = 1.4;
Pi = Pc*1.18; % manifold pressure, psi
R_univ = 10.73159; %psi*ft3/lbmol*R
mdot_ox = CdA_lox_inj/12*sqrt(2*g*rho_ox*(Pi-Pc));
mdot_fuel = CdA_fuel_inj/12*sqrt(2*g*rho_fuel*(Pi-Pc));

%create scattered interpolants
F_T = scatteredInterpolant(CEA.p,CEA.OF,CEA.t);
F_mw = scatteredInterpolant(CEA.p,CEA.OF,CEA.mw);
F_gam = scatteredInterpolant(CEA.p,CEA.OF,CEA.gam);


for i = 1:length(t)
    if i > 1
        Pi1 = (12*mdot_ox/CdA_lox_inj)^2*1/(2*g*rho_ox)+Pc;
        Pi2 = (12*mdot_fuel/CdA_fuel_inj)^2*1/(2*g*rho_fuel)+Pc;
        Pi = (Pi1+Pi2)/2;
    end
    R = R_univ/F_mw(Pc,OF);
    gam = F_gam(Pc,OF);
    mdot_ex = Pc*At/12*sqrt(g*gam/(R*T_c))*(1+(gam-1)/2)^(-(gam+1)/(2*(gam-1)));
    dp_c = (R*T_c/Vc*(mdot_ox+mdot_fuel-mdot_ex))*dt;
    Pc = Pc + dp_c;
    if Pc < 0 || ~isreal(Pc)
        Pc = 0;
    end
    
    % recalculate OF w new PC
    mdot_ox = CdA_lox_inj/12*sqrt(2*g*rho_ox*(Pi-Pc));
    mdot_fuel = CdA_fuel_inj/12*sqrt(2*g*rho_fuel*(Pi-Pc));
    OF = mdot_ox/mdot_fuel;
    if real(OF) < 0.6
        OF = 0.6;
    elseif real(OF) > 1.7
        OF = 1.7;
    elseif ~isreal(OF)
        OF = 1.4;
    end
    T_c = F_T(Pc,OF)*1.8;
    
    Vdot = mdot_ox/rho_ox + mdot_fuel/rho_fuel;
    mdot_press = Q_tot*Pt/(R_HE*T_press);
    
    Pc_arr(i) = Pc;
    T_arr(i) = T_c;   
    OF_arr(i) = OF;
    
end
%%
figure
plot(t,Pc_arr)
ylabel('Chamber Pressure psi','LineWidth',2)
xlabel('Time s')
xlim([0 0.03])
grid on

%%
figure
plot(t,T_arr./1.8)
ylabel('Temperature K','LineWidth',2)
xlabel('Time s')
xlim([0 0.03])
grid on

%%
figure
plot(t,OF_arr)
ylabel('Temperature K','LineWidth',2)
xlabel('Time s')
xlim([0 0.03])
grid on