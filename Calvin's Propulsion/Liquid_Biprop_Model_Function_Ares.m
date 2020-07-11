%function [t,FT,I,Isp,mdot_tot,OF,p_p,p_f,p_o,p_c,m_f,m_o,p_fmax,p_omax,p_cmax,unstabI,unstabG,i] = Liquid_Biprop_Model_Function(rho_f,rho_o,OF_ref,c_star_ref,V_p,V_f,V_o,V_fpi,V_opi,p_pi,p_fset,p_oset,Cd_f,Cd_o,A_f,A_o,A_t)
close all; clear all; clc;
addpath('Combustion Data')
%% Regulator Set Pressure
p_fset = 500;               % fuel regulator initial set pressure [psia]
p_oset = 500;               % oxidizer regulator initial set pressure [psia]

%% Tank Volumes
V_f = 0.360;                % fuel tank volume [ft^3]
V_o = 0.422;                % oxidizer tank volume [ft^3]
V_p = 124/12^3;

%% Initial Propellant Volumes
V_fpi = V_f-0.001;
V_opi = V_o-0.001;

%% Initial Tank Pressures
p_pi = 4500;               % initial pressurant tank pressure [psia]
p_fi = p_fset;             % initial fuel tank pressure [psia]
p_oi = p_oset;             % initial oxidizer tank pressure [psia]

%% Initial Tank Temperatures
T_pi = 528;                 % initial pressurant tank temperature [degR]
T_fui = 528;                % initial fuel tank temperature [degR]
T_oui = 528;                % initial oxygen tank temperature [degR]

%% Discharge Coefficient, CdA (Liquid/Injector)
Cd_f = 0.5;                 % fuel discharge coefficient [dim]
Cd_o = 0.5;                 % oxidizer discharge coefficient [dim]
A_f = 0.027;                % fuel flow area [in^2]
A_o = 0.050;                % oxidizer flow area [in^2]

%% Flow Coefficient, Cv (Gas/Regulator)
p_hyst = 5;
Cv_Rf = 0.02;                % flow coefficient of fuel regulator [dim]
Cv_Ro = 0.02;                % flow coefficient of oxidizer regulator [dim]
SPE_Rf = 1/100;             % supply pressure effect of fuel reg [psi increase outlet/psi decrease inlet]
SPE_Ro = 1/100;             % supply pressure effect of ox reg [psi increase outlet/psi decrease inlet]

%% Engine Performance Properties
cstar_eff = 90;             % combustion (c*) efficiency [%]
cT = 1.4;                   % thrust coefficient, nozzle performance [dim]
percentdp_crit = 20;        % critical pressure drop across injector [%]
A_t = (1.5/2)^2*pi;         % throat area [in^2]

%% Propellant Properties
load('RP1LOx.mat');
rho_f = rho_fu_raw/32.2;    % fuel density [slug/ft^3]
rho_o = rho_ox_raw/32.2;    % oxidizer density [slug/ft^3]

%% Helium Properties
gam_p = 1.66;               % specific heat ratio for He [dim]
R_p = 12421;                % specific gas constant for He [ft*lbf/slug*degR]
G_p = 0.138;                % specific gravity for He [dim]
c_v = R_p/(gam_p - 1);      % specific heat capacity under constant volume [ft*lbf/slug*degR]
c_p = R_p*gam_p/(gam_p - 1);% specific heat capacity under constant pressure [ft*lbf/slug*degR]

a_p_SI = 3.46E-3;            % van der Waal's coefficient a for He [Pa m^6]
a_p = a_p_SI*(0.020885)*(3.2808)^6;  % vdW coefficient a for He [lbf/ft^2 ft^6]
b_p_SI = 23.71E-6;           % van der Waal's coefficient b for He [m^3/mol]
b_p = b_p_SI*(3.2808)^3;     % van der Waal's coefficient b for He [ft^3/mol]

%% Ulage Volume Calc
V_fui = V_f-V_fpi;           % initial fuel ullage volume [ft^3]
V_oui = V_o-V_opi;           % initial oxidizer ullage volume [ft^3]
m_fi = rho_f*V_fpi;          % initial fuel mass [slug]
m_oi = rho_o*V_opi;          % initial oxidizer mass [slug]

%% Regulator Flow Rate Calc
Cvcc = 22.67/60;             % Cv correction constant for std ft^3/s

%% Simulation Parameters
t_final = 60;
dt = 0.001;
N_sim = ceil(t_final/dt);

%% Array Prep
t = zeros(1,N_sim);
p_p = zeros(1,N_sim);
p_f = zeros(1,N_sim);
p_o = zeros(1,N_sim);
p_c = zeros(1,N_sim);
FT = zeros(1,N_sim);
T_p = zeros(1,N_sim);
T_fu = zeros(1,N_sim);
T_ou = zeros(1,N_sim);
V_fu = zeros(1,N_sim);
V_ou = zeros(1,N_sim);
rho_p = zeros(1,N_sim);
rho_fu = zeros(1,N_sim);
rho_ou = zeros(1,N_sim);
m_p = zeros(1,N_sim);
m_fu = zeros(1,N_sim);
m_ou = zeros(1,N_sim);
m_ptot = zeros(1,N_sim);
m_f = zeros(1,N_sim);
m_o = zeros(1,N_sim);
dp_Rf = zeros(1,N_sim);
dp_Ro = zeros(1,N_sim);
mdot_Rf = zeros(1,N_sim);
mdot_Ro = zeros(1,N_sim);
mdot_pR = zeros(1,N_sim);
mdot_f = zeros(1,N_sim);
mdot_o = zeros(1,N_sim);
mdot_tot = zeros(1,N_sim);
%dp_f = zeros(1,N_sim);
%dp_o = zeros(1,N_sim);
OF = zeros(1,N_sim);
%Isp = zeros(1,N_sim);
cstar = zeros(1,N_sim);
cstar_ideal = zeros(1,N_sim);

for j=1:N_sim
   A{j} = zeros(14);
   MDOT{j} = zeros(14,1);
   X{j} = zeros(14,1);
   dX{j} = zeros(14,1);
   %PcPolyroots{j} = zeros(1,3);
end

t(1) = 0;
p_p(1) = p_pi;
p_f(1) = p_fi;
p_o(1) = p_oi;
p_c(1) = 15;                             % chamber pressure [psia]
FT(1) = 0;
T_p(1) = T_pi;
T_fu(1) = T_fui;
T_ou(1) = T_oui;
V_fu(1) = V_fui;                        % fuel ullage volume [ft^3]
V_ou(1) = V_oui;                        % oxidizer ullage volume [ft^3]
rho_p(1) = p_pi/(R_p*T_pi)*144;         % pressurant density [slug/ft^3]
rho_pwdW(1) = RhoVDW(p_pi*144,T_pi,R_p,a_p,b_p);
rho_fu(1) = p_fi/(R_p*T_fui)*144;       % fuel ullage density [slug/ft^3]
rho_ou(1) = p_oi/(R_p*T_oui)*144;       % oxidizer ullage density [slug/ft^3]
m_p(1) = rho_p(1)*V_p;                  % pressurant mass [slug]
m_fu(1) = rho_fu(1)*V_fu(1);            % fuel ullage mass [slug]
m_ou(1) = rho_ou(1)*V_ou(1);            % oxidizer ullage mass [slug]
m_f(1) = m_fi;
m_o(1) = m_oi;

unstabI = 0;
unstabG = 0;
i = 1;

while (m_f(i) > 0) && (m_o(i) > 0) && t(i) < t_final
    
    p_fseteff = p_fset + SPE_Rf*(p_pi - p_p(i));                    % fuel effective regulator set pressure [psia]
    dp_Rf(i) = p_p(i) - p_f(i);                                     % fuel regulator pressure differential [psia]
    if p_f(i) >= p_fseteff                                          % regulator shut off
        mdot_Rf(i) = 0;
    elseif p_p(i)/p_f(i) < 2.05                                     % unchoked flow
        StdVdot_Rf = Cvcc*Cv_Rf*p_p(i)*(1 - 2*dp_Rf(i)/(3*p_p(i)))...
        *sqrt(dp_Rf(i)/(p_p(i)*G_p*T_p(i)));                        % fuel pressurant volume flowrate [std ft^3/s]
        mdot_Rf(i) = StdVdot_Rf*14.7*144/(R_p*528);                 % fuel pressurant mass flowrate [slug/s]
    elseif p_p(i)/p_f(i) >= 2.05                                    % choked flow
        StdVdot_Rf = 0.471*Cvcc*Cv_Rf*p_p(i)*sqrt(1/(G_p*T_p(i)));  % fuel pressurant volume flowrate [std ft^3/s]
        mdot_Rf(i) = StdVdot_Rf*14.7*144/(R_p*528);                 % fuel pressurant mass flowrate [slug/s]
    else
        error('pressure error');
    end


%     if i == 1
%         stdq_Rf = RegGasFlow(Cv_Rf,0,p_p(i),p_f(i),p_fset,p_hyst,G_p,T_p(i),gam_p);
%         stdq_Ro = RegGasFlow(Cv_Ro,0,p_p(i),p_o(i),p_oset,p_hyst,G_p,T_p(i),gam_p);
%     else
%         stdq_Rf = RegGasFlow(Cv_Rf,stdq_Rf,p_p(i),p_f(i),p_fset,p_hyst,G_p,T_p(i),gam_p);
%         stdq_Ro = RegGasFlow(Cv_Ro,stdq_Ro,p_p(i),p_o(i),p_oset,p_hyst,G_p,T_p(i),gam_p);
%     end
%     mdot_Rf(i) = stdq_Rf*14.7*144/(R_p*528)/60;
%     mdot_Ro(i) = stdq_Ro*14.7*144/(R_p*528)/60;
    
    p_oseteff = p_oset + SPE_Ro*(p_pi - p_p(i));                    % oxidizer effective regulator set pressure [psia]
    dp_Ro(i) = p_p(i) - p_f(i);                                     % oxidizer regulator pressure differential [psia]
    if p_o(i) >= p_oseteff                                          % regulator shut off
        mdot_Ro(i) = 0;
    elseif p_p(i)/p_o(i) < 2.05                                     % unchoked flow
        StdVdot_Ro = Cvcc*Cv_Ro*p_p(i)*(1 - 2*dp_Ro(i)/(3*p_p(i)))...
        *sqrt(dp_Ro(i)/(p_p(i)*G_p*T_p(i)));                        % oxidizer pressurant volume flowrate [std ft^3/s]
        mdot_Ro(i) = StdVdot_Ro*14.7*144/(R_p*528);                 % oxidizer pressurant mass flowrate [slug/s]    
    elseif p_p(i)/p_o(i) >= 2.05                                    % choked flow
        StdVdot_Ro = 0.471*Cvcc*Cv_Ro*p_p(i)*sqrt(1/(G_p*T_p(i)));  % oxidizer pressurant volume flowrate [std ft^3/s]
        mdot_Ro(i) = StdVdot_Ro*14.7*144/(R_p*528);                 % oxidizer pressurant mass flowrate [slug/s]
    else
        error('pressure error');
    end
    
    mdot_pR(i) = mdot_Rf(i) + mdot_Ro(i);
    
    if i == 1
        OF(i) = (Cd_o*A_o)/(Cd_f*A_f)*sqrt(rho_o/rho_f);
    else
        OF(i) = mdot_o(i-1)/mdot_f(i-1);
    end
    
    cstar_ideal(i) = interp1(OF_ref,c_star_ref,OF(i),'spline');          % ideal characteristic velocity [ft/s]
    cstar(i) = cstar_ideal(i)*cstar_eff*.01;                       % real characteristic velocity [ft/s]
    
    p_cinf = linspace(0,p_f(i),100);
    p_cino = linspace(0,p_o(i),100);
    dP_f = p_f(i) - p_cinf;
    dP_o = p_o(i) - p_cino;
    mdot_fM = Cd_f*A_f*sqrt(2*rho_f*dP_f/144);
    mdot_oM = Cd_o*A_o*sqrt(2*rho_o*dP_o/144);
    mdot_totM = mdot_fM + mdot_oM;
    p_cout = mdot_totM*cstar(i)/A_t;
    p_c(i) = interp1(p_cout-p_cinf,p_cinf,0);
    mdot_f(i) = interp1(p_cout-p_cinf,mdot_fM,0);
    mdot_o(i) = interp1(p_cout-p_cinf,mdot_oM,0);
    mdot_tot(i) = mdot_f(i) + mdot_o(i);
    m_f(i+1) = m_f(i) - mdot_f(i)*dt;                           % iterate fuel propellant mass in tank [slug]
    m_o(i+1) = m_o(i) - mdot_o(i)*dt;                           % iterate oxidizer propellant mass in tank [slug]
    
    dE_p(i) = c_p*(-mdot_pR(i)*T_p(i));
    dE_fu(i) = c_p*(mdot_Rf(i)*T_p(i));
    dE_ou(i) = c_p*(mdot_Ro(i)*T_p(i));
    
%     if i == 1
%         OF(i) = (Cd_o*A_o)/(Cd_f*A_f)*sqrt(rho_o/rho_f);
%     else
%         OF(i) = mdot_o(i-1)/mdot_f(i-1);
%     end
%     
%     cstar_ideal(i) = interp1(OF_ref,c_star_ref,OF(i));          % ideal characteristic velocity [ft/s]
%     cstar(i) = cstar_ideal(i)*cstar_eff*.01;                       % real characteristic velocity [ft/s]
%     PcPolyroots{i} = roots([A_t^2/cstar(i)^2,...                   % roots of chamber pressure polynomial [psia]
%          (2*rho_f*(Cd_f*A_f)^2 + 2*rho_o*(Cd_o*A_o)^2)/144,...
%          (-p_f(i)*(2*rho_f*(Cd_f*A_f)^2)-p_o(i)*(2*rho_o*(Cd_o*A_o)^2))/144]);
%      p_c2(i) = max(PcPolyroots{i}(PcPolyroots{i} >= 0));          % chamber pressure [psia]
%     mdot_veri(i) = p_c(i)*A_t/cstar(i);
%     
%     dp_f(i) = p_f(i) - p_c(i);                                  % fuel injector pressure differential [psia]
%     mdot_f(i) = (Cd_f*A_f)*sqrt(2*rho_f*dp_f(i)/144);           % fuel propellant mass flowrate [slugs/s]
%     m_f(i+1) = m_f(i) - mdot_f(i)*dt;                           % iterate fuel propellant mass in tank [slug]
%     dp_o(i) = p_o(i) - p_c(i);                                  % oxidizer injector pressure differential [psia]
%     mdot_o(i) = (Cd_o*A_o)*sqrt(2*rho_o*dp_o(i)/144);           % oxidizer propellant mass flowrate [slugs/s]
%     m_o(i+1) = m_o(i) - mdot_o(i)*dt;                           % iterate oxidizer propellant mass in tank [slug]
%     mdot_tot(i) = mdot_f(i) + mdot_o(i);
%     
    FT(i) = p_c(i)*A_t*cT;
    FT_check = mdot_tot(i)*cstar(i)*cT;
    if abs(FT(i) - FT_check) > 1
        fprintf('Thrust values don''t match!')
    end
    %Isp(i) = FT(i)/(mdot_tot(i)*32.2);

    
    A{i} = [0,0,-V_p,1,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,0,-V_fu(i),1,-rho_fu(i),0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,0,-V_ou(i),1,-rho_ou(i);
        0,0,0,-1,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,1,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,0,0,1,0;
        0,0,0,0,0,0,0,0,rho_f,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,0,0,0,rho_o;
        1,-R_p*rho_p(i),-R_p*T_p(i),0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,1,-R_p*rho_fu(i),-R_p*T_fu(i),0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,1,-R_p*rho_ou(i),-R_p*T_ou(i),0,0;
        0,-m_p(i)/(gam_p*T_p(i)),0,-1/gam_p,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,m_fu(i)/(gam_p*T_p(i)),0,T_fu(i)/(gam_p*T_p(i)),(gam_p-1)*p_f(i)*144/(gam_p*R_p*T_p(i)),0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,m_ou(i)/(gam_p*T_p(i)),0,T_ou(i)/(gam_p*T_p(i)),(gam_p-1)*p_o(i)*144/(gam_p*R_p*T_p(i))];
    
    MDOT{i} = [0;0;0;mdot_pR(i);mdot_Rf(i);mdot_Ro(i);mdot_f(i);mdot_o(i);0;0;0;mdot_pR(i);mdot_Rf(i);mdot_Ro(i)];
    
    dX{i} = A{i}\MDOT{i};
    X{i} = [p_p(i)*144;T_p(i);rho_p(i);m_p(i);p_f(i)*144;T_fu(i);rho_fu(i);m_fu(i);V_fu(i);p_o(i)*144;T_ou(i);rho_ou(i);m_ou(i);V_ou(i)];
    X{i+1} = X{i} + dX{i}*dt;
    
    
    if ((p_f(i)-p_c(i)) <= p_c(i)*percentdp_crit*0.01) || ((p_o(i)-p_c(i)) <= p_c(i)*percentdp_crit*0.01)
        unstabI = 1;
        fprintf('Unstable Across Injector\n')
        break
    elseif m_f(i) < 0
        fprintf('Out of Fuel\n')
    elseif m_o(i) < 0
        fprintf('Out of Oxidizer\n')
    elseif i == N_sim - 1
        fprintf('Burn Time Too Long\n')
        unstabG = 1;
    end
    
    i = i + 1;

    [p_p(i),T_p(i),rho_p(i),m_p(i),p_f(i),T_fu(i),rho_fu(i),m_fu(i),V_fu(i),p_o(i),T_ou(i),rho_ou(i),m_ou(i),V_ou(i)] = matsplit(X{i});
    m_ptot(i) = m_p(i) + m_fu(i) + m_ou(i);
    p_p(i) = p_p(i)/144;
    p_f(i) = p_f(i)/144;
    p_o(i) = p_o(i)/144;
    
    t(i) = t(i-1) + dt;                 % iterate time [s]
end

if i > 1
    i = i - 1;
end

H_p = m_p.*c_p.*T_p;
H_fu = m_fu.*c_p.*T_fu;
H_ou = m_ou.*c_p.*T_ou;
H_tot = H_p+H_fu+H_ou;

U_p = m_p.*c_v.*T_p;
U_fu = m_fu.*c_v.*T_fu;
U_ou = m_ou.*c_v.*T_ou;
U_tot = U_p+U_fu+U_ou;


I = cumtrapz(FT)*dt;
Isp = I(i)/(m_f(1)+m_o(1))/32.2;

p_fmax = max(p_f(1:i));
p_omax = max(p_o(1:i));
p_cmax = max(p_c(1:i));

fprintf('Total Impulse: %.2f lbf-s\n',I(i))
fprintf('Average Thrust: %.2f lbf\n',mean(FT(1:i)))
fprintf('Isp: %.2f s\n',Isp)
fprintf('Burn Time: %.2f s\n',t(i))
fprintf('\n')

Lw = 3;
Fs = 14;

figure
set(gca,'fontsize',Fs)
plot(t(1:i),FT(1:i),'Linewidth',Lw)
xlabel('Time [s]')
ylabel('Thrust [lbf]')
title('Thrust Curve')

figure
set(gca,'fontsize',Fs)
plot(t(1:i),OF(1:i),'Linewidth',Lw)
xlabel('Time [s]')
ylabel('OF [dim]')

% figure
% set(gca,'fontsize',Fs)
% plot(t(1:i),Isp(1:i),'Linewidth',Lw)
% xlabel('Time [s]')
% ylabel('Specific Impulse, Isp [s]')

figure
set(gca,'fontsize',Fs)
hold on
plot(t(1:i),p_p(1:i),'Linewidth',Lw)
plot(t(1:i),p_f(1:i),'Linewidth',Lw)
plot(t(1:i),p_o(1:i),'Linewidth',Lw)
plot(t(1:i),p_c(1:i),'Linewidth',Lw)
plot(t(1:i),p_c(1:i)*0.8,'Linewidth',Lw)
ylim([0,6000])
hold off
xlabel('Time [s]')
ylabel('Pressure [psia]')
legend('Pressurant','Fuel','Oxidizer','Chamber','Critical')

figure
set(gca,'fontsize',Fs)
hold on
plot(t(1:i),m_p(1:i),'Linewidth',Lw)
plot(t(1:i),m_fu(1:i),'Linewidth',Lw)
plot(t(1:i),m_ou(1:i),'Linewidth',Lw)
plot(t(2:i),m_ptot(2:i),'Linewidth',Lw)
hold off
xlabel('Time [s]')
ylabel('Gas Mass [slug]')
legend('Pressurant','Fuel','Oxidizer','Total')

figure
set(gca,'fontsize',Fs)
hold on
plot(t(1:i),dE_p(1:i),'Linewidth',Lw)
plot(t(1:i),dE_fu(1:i),'Linewidth',Lw)
plot(t(1:i),dE_ou(1:i),'Linewidth',Lw)
plot(t(1:i),dE_p(1:i)+dE_fu(1:i)+dE_ou(1:i),'Linewidth',Lw)
hold off
xlabel('Time [s]')
ylabel('Energy [slug]')
legend('Pressurant','Fuel','Oxidizer','Total')

figure
set(gca,'fontsize',Fs)
hold on
plot(t(1:i),m_f(1:i),'Linewidth',Lw)
plot(t(1:i),m_o(1:i),'Linewidth',Lw)
hold off
xlabel('Time [s]')
ylabel('Propellant Mass [slug]')
legend('Fuel','Oxidizer')

figure
set(gca,'fontsize',Fs)
hold on
plot(t(1:i),mdot_pR(1:i),'Linewidth',Lw)
plot(t(1:i),mdot_Rf(1:i),'Linewidth',Lw)
plot(t(1:i),mdot_Ro(1:i),'Linewidth',Lw)
hold off
xlabel('Time [s]')
ylabel('Gas Mass Flow Rate [slug/s]')
legend('Pressurant','Fuel','Oxidizer')

figure
set(gca,'fontsize',Fs)
hold on
plot(t(1:i),T_p(1:i),'Linewidth',Lw)
plot(t(1:i),T_fu(1:i),'Linewidth',Lw)
plot(t(1:i),T_ou(1:i),'Linewidth',Lw)
hold off
xlabel('Time [s]')
ylabel('Temperature [psia]')
legend('Pressurant','Fuel','Oxidizer')

function rho = RhoVDW(p, T, R, a, b)
%PVDW VanDerWaal's equation of state for density, gas values
    vdw_roots = roots([-a*b, a, - p*b - R*T, p]);
    real_vdw_roots = [];
    for ii = 1:length(vdw_roots)
        if isreal(vdw_roots(ii))
            real_vdw_roots(end+1) = vdw_roots(ii);
        end
    end
    rho = min(real_vdw_roots);
end
