close all;
%% Import and Scale/Offset Data

dat = lvm_import('static_fire_2_22.lvm',1); % import labview data file
dat = dat.Segment1.data; % get actual data

time = dat(:,1);
oxman = 242.47*dat(:,3) - 237.6;
fuelman = dat(:,13);
% CCRaw =  245.23*dat(:,5)- 240.52;
CCRaw =  245.23*dat(:,5) - 240.52;
oxtank = 252.39*dat(:,6) - 252.92;
press = dat(:,11);
% pnuem = 61.045*PneumaticsRaw - 63.905;
fueltank = oxtank; % fuel tank not aquired for this fire
Thrust = dat(:,9);
Thrust = Thrust*0.9191-23.1582;

ind = 8301:8824; % data index range for burn

% CC = CCRaw;

% Zero the data set time
t0 = 1905000;
time0 = (time - t0)/1000;
burnind = (ind(1)+30):(ind(end)-200); % for processing data relevent to the burn;

%% Engine Properties
At = 1.77; % throat area, in2
Ae = pi/4*2.81^2; % exit area, in2
C_tau = 1.4;
g = 32.174;
rho_lox = 72.2; % lb / ft3
rho_fuel = 53.4; % lb / ft3

CdA_lox_inj = 0.0453; % in ** 2
CdA_fuel_inj = 0.0339; % in ** 2
CdA_lox_lines = 0.0410;
CdA_fuel_lines = 0.0564;
CdA_lox = 1/(sqrt(1/CdA_lox_inj^2+1/CdA_lox_lines^2));
CdA_fuel = 1/(sqrt(1/CdA_fuel_lines^2+1/CdA_fuel_inj^2));

derivedPC = Thrust(ind)./(At.*1.4);

%% Predict C_tau theoeretical
load CEADAT_LOXETH_93WW.mat
F_T = scatteredInterpolant(CEA.p,CEA.OF,CEA.t);
F_mw = scatteredInterpolant(CEA.p,CEA.OF,CEA.mw);
F_gam = scatteredInterpolant(CEA.p,CEA.OF,CEA.gam);
F_cstar = scatteredInterpolant(CEA.p,CEA.OF,CEA.cstar);

gam = F_gam(250, 1.4); % get specific heat ratio from interpolated CEA
% Solve for exit mach number from area ratio
M = 3; % guess mach number
err = 1;
c1 = (gam + 1)/2; c2 = gam - 1;

while abs(err) > 0.001
    Arat = c1^(-c1/c2)*((1+c2/2*M^2)^(c1/c2))/M;
    err = Arat - Ae/At;
    M = M-err*0.1;
end 

% Solve for exit pressure based on mach number
Pe = (14.7+CCRaw).*(1+c2/2*M^2)^(-gam/c2); % exit pressure,psi
Pa = 14.7; % ambient pressure
C_tau_theo = sqrt(2*gam^2/c2*(1/c1)^(2*c1/c2)*(1-(Pe./(14.7+CCRaw)).^(c2/gam)))+(Pe-Pa).*Ae./((14.7+CCRaw).*At);

% Get rid of nonsense startup and ending transient values
C_tau_theo = [ones(ind(1)+25-1,1)*C_tau_theo(ind(1)+25);...
              C_tau_theo((ind(1)+25):(ind(end)-100));                        ...
              ones(length(Thrust)-(ind(end)-100),1)*C_tau_theo(ind(end)-100)];
% CC = Thrust ./ (At.*C_tau_theo); % PC derived from thrust and thrust coefficient
CC = Thrust ./ (At.*1.41); % PC derived from thrust and thrust coefficient

%% Thrust Coefficient
c_tau = Thrust ./ (CC .* At);
figure
plot(time0(ind),C_tau_theo(ind),'LineWidth',2)
title('C\tau Theoretical')
xlabel('Time (s)')
ylabel('Thrust Coefficient, %')
legend('C\tau');
ylim([1.2,1.6])
xlim([0 20]);
set(gca, 'FontSize', 20, 'FontWeight', 'bold')
grid on

%% Calculate Mass Flow Rates
mdot_lox_1 = real(CdA_lox/12 * sqrt(2.*g.*rho_lox.*(oxtank-CC)));
mdot_lox_2 = real(CdA_lox_inj/12 * sqrt(2.*g.*rho_lox.*(oxman-CC)));
mdot_lox_3 = real(CdA_lox_lines/12 * sqrt(2.*g.*rho_lox.*(oxtank-oxman)));
mdot_fuel_1 = real(CdA_fuel/12 * sqrt(2.*g.*rho_fuel.*(fueltank-CC)));
mdot_fuel_2 = real(CdA_fuel_inj/12 * sqrt(2.*g.*rho_fuel.*(fuelman-CC)));
mdot_fuel_3 = real(CdA_fuel_lines/12 * sqrt(2.*g.*rho_fuel.*(fueltank-fuelman)));

% Compare Mass Flow Calculations
figure
subplot(2,1,1)
plot(time0(ind), mdot_lox_1(ind), time0(ind), mdot_lox_2(ind), time0(ind), mdot_lox_3(ind))
legend('Ox: Tank-CC','Ox: Manifold-CC','Ox: Tank-Man')
ylim([0 3])
subplot(2,1,2)
plot(time0(ind), mdot_fuel_1(ind),time0(ind), mdot_fuel_2(ind),time0(ind), mdot_fuel_3(ind))
ylim([0 2])
xlabel('Time (s)')
ylabel('Mass Flow Rate (lbm/s)')
legend('Fuel: Tank-CC','Fuel: Manifold-CC','Fuel: Tank-Man')

%% Plot best mass flow rates
figure
plot(time0(ind),smoothdata(mdot_lox_1(ind),'gaussian',8), time0(ind), smoothdata(mdot_fuel_2(ind),'gaussian',8),'LineWidth',2)
legend('Ox','Fuel')
ylim([0 3.5])
xlabel('Time (s)')
ylabel('Mass Flow Rate (lbm/s)')
grid on

%% Chamber Pressure Comparison
figure
plot(time0(ind), CC(ind), time0(ind), CCRaw(ind), 'LineWidth', 2);
xlim([0 20]);
title('Chamber Pressure Calculated vs. Measured', 'FontSize', 16)
ylabel('Pressure (psig)')
legend('Calculated from Thrust', 'Measured');
xlabel('Time (s)')
set(gca, 'FontSize', 17, 'FontWeight', 'bold')

%% Tank Pressures
figure
subplot(2,1,1);
plot(time0(ind),press(ind),'LineWidth',2)
grid on
title('Tank Pressures vs Time')
ylabel('Pressure (psig)')
legend('Pressurant Tank');
set(gca, 'FontSize', 18, 'FontWeight', 'bold')
subplot(2,1,2);
plot(time0(ind),oxtank(ind),time0(ind),fueltank(ind),'LineWidth',2);
ylabel('Pressure (psig)')
xlabel('Time (s)')
legend('Oxidizer Tank (psig)','Fuel Tank (psig)')
set(gca, 'FontSize', 18, 'FontWeight', 'bold')
grid on
%% Tank Pressures
figure
subplot(2,1,1);
plot(time0,press,'LineWidth',2)
grid on
title('Tank Pressures vs Time')
ylabel('Pressure (psig)')
legend('Pressurant Tank');
xlim([-150, 30])
set(gca, 'FontSize', 18, 'FontWeight', 'bold')
subplot(2,1,2);
plot(time0,oxtank,time0,fueltank,'LineWidth',2);
ylabel('Pressure (psig)')
xlabel('Time (s)')
legend('Oxidizer Tank (psig)','Fuel Tank (psig)')
set(gca, 'FontSize', 18, 'FontWeight', 'bold')
xlim([-150, 30])
grid on

%% All Presures
figure
plot(time0(ind),oxtank(ind),time0(ind),fueltank(ind),time0(ind),oxman(ind),time0(ind),fuelman(ind), ...
    time0(ind), CC(ind),time0(ind), CCRaw(ind), 'LineWidth',2);
ylabel('Pressure (psig)')
xlabel('Time (s)')
legend('Oxidizer Tank','Fuel Tank','Ox Manifold','Fuel Manifold','Derived PC','Raw PC')
set(gca, 'FontSize', 18, 'FontWeight', 'bold')
grid on
xlim([-1 17])
%% Shutdown Transient
figure
hold on
plot(time0(ind),oxman(ind),'g',time0(ind),fuelman(ind),'b', ...
    time0(ind), CCRaw(ind),'m', 'LineWidth',2);
ylabel('Pressure (psig)')
xlabel('Time (s)')
set(gca, 'FontSize', 18, 'FontWeight', 'bold')
ylim([-3,300])
yyaxis right
plot(time0(ind), Thrust(ind),'LineWidth',2,'Color','r')
ylim([0,50])
xlim([12 24])
set(gca, 'FontSize', 18, 'FontWeight', 'bold')
legend('Ox Manifold','Fuel Manifold','Chamber Pressure','Thrust')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
grid on

%% Thrust and Chamber Pressure
figure
hold on
plot(time0(ind),CCRaw(ind),'LineWidth',2)
ylabel('Pressure (psig)')
xlabel('Time(s)')
xlim([0 17])
ylim([0 300])
set(gca, 'FontSize', 20, 'FontWeight', 'bold')
yyaxis right
plot(time0(ind),Thrust(ind),'LineWidth',2);
title('Thrust and Chamber Pressure')
ylabel('Thrust (lbf)')
xlabel('Time(s)')
xlim([0 17])
ylim([0 700])
set(gca, 'FontSize', 20, 'FontWeight', 'bold')
legend('Chamber Pressure','Thrust')
grid on
%% Manifold Pressures
figure
plot(time0(ind),oxman(ind),time0(ind),fuelman(ind),'LineWidth',2)
title('Manifold Pressures')
xlabel('Time (s)')
ylabel('Pressure (psig)')
xlim([0 20]);
ylim([0 400]);
legend('Ox','Fuel Manifold');
set(gca, 'FontSize', 18, 'FontWeight', 'bold')
grid on

%% Injector Stiffness w.r.t Time
oxdppc = (oxman-CC)./CC*100;
fueldppc = (fuelman-CC)./CC*100;
figure
hold on
title('Injector Pressure Drop Normalized by Chamber Pressure')
plot(time0(ind),smoothdata(oxdppc(ind),'gaussian',6),time0(ind),smoothdata(fueldppc(ind),'gaussian',6),'LineWidth',2)
instabilitylimit = 13;
yline(instabilitylimit,'--','LineWidth',2);
s = patch([0,0,time0(end),time0(end)],[0,instabilitylimit,instabilitylimit,0],[0.2,0,0]+0.8);
s.FaceVertexAlphaData = 0.001;
s.FaceAlpha = 'flat';
legend('Ox DP*','Fuel DP*')
ylabel('Normalized Injector Pressure Drop Percentage (%)')
xlabel('Time (s)')
set(gca, 'FontSize', 20, 'FontWeight', 'bold')
xlim([5,13]);
ylim([0 70]);
grid on
%% Mixture Ratio
OF = abs(mdot_lox_1)./abs(mdot_fuel_1); % Mixture ratio,Ox: Tank-CC Fuel: Man-CC
OF_optimal = 1.45;

figure
plot(time0(burnind),smoothdata(OF(burnind),'gaussian',8),'LineWidth',2)
hold on
yline(OF_optimal,'--','LineWidth',2);
title('Mixture Ratio')
xlabel('Time (s)')
ylabel('Mixture Ratio')
legend('Measured','Ideal');
ylim([0.7,1.8]);
xlim( [-1 15]);
grid on
set(gca, 'FontSize', 20, 'FontWeight', 'bold')

%% Combustion Efficiency

cstar_ideal = F_cstar(220, 1.3)*3.28; % ft / s
mdot = abs(mdot_fuel_1)+abs(mdot_lox_1);
cstar = CC.* At .* g ./ mdot;
figure
hold on
plot(time0(ind),smoothdata(cstar(ind),'gaussian',5),'LineWidth',2)
yline(cstar_ideal,'--','LineWidth',2);
title('C^*')
xlabel('Time (s)')
ylabel('C^*, ft / s')
legend('Measured','Ideal');
xlim([0 20]);
set(gca, 'FontSize', 15, 'FontWeight', 'bold')
grid on


cstar_eff = 100.*cstar./cstar_ideal;
figure
plot(time0(ind),smoothdata(cstar_eff(ind),'gaussian',5),'LineWidth',2)
title('C^* Efficiency')
xlabel('Time (s)')
ylabel('C^* Efficiency, %')
xlim([0 20]);
ylim([0,105]);
set(gca, 'FontSize', 15, 'FontWeight', 'bold')
grid on

%% Adiabatic Flame Temperature
Mw = F_mw(250,1.3); % kg / kmol
R = 8314;
n1 = gam + 1;
n2 = gam - 1;
T0 = (cstar/3.28).^2.*gam*(2/n1)^(n1/n2).*Mw/R;
T0_ideal = F_T(250,1.4); % K
T0_95 = (cstar_ideal/3.28*0.75).^2.*1.1346*(2/n1)^(n1/n2).*Mw/R;
figure
hold on
plot(time0(ind),smoothdata(T0(ind),'gaussian',5),'LineWidth', 2)
yline(T0_ideal,'--','LineWidth', 2)
yline(T0_95,'--','LineWidth', 2,'Color','r')
xlabel('Time (s)')
ylabel('Adiabatic Flame Temperature, K')
legend('Measured','Ideal 95% Ethanol','95% Ethanol, 75% C^*_{eff}');
grid on;
xlim([-1,15])
ylim([500,3500])
set(gca, 'FontSize', 20, 'FontWeight', 'bold')
%% Specific Impulse
Isp = cstar.*c_tau./g;
Isp_ideal = cstar_ideal*1.4/g;
figure
hold on
plot(time0(ind),smoothdata(Isp(ind),'gaussian',5),'LineWidth',2)
yline(Isp_ideal,'--','LineWidth',2);
title('Specific Impulse')
xlabel('Time (s)')
ylabel('I_s_p, s')
xlim([0 18]);
grid on
set(gca, 'FontSize', 15, 'FontWeight', 'bold')
%% Supply Pressure Effect Curve
SPE_ind = 5732:5996; 
x = press(SPE_ind);
y = oxtank(SPE_ind);
X = [ones(length(x),1) x];
b = X\y;
SPE = -b(2); % Slope of line, supply pressure effect

figure
hold on
plot(x,y)
plot(x,b(2).*x+b(1),'b--','LineWidth',2)
xlim([1500 4000])
ylim([350 400])
xl = xlim;
yl = ylim;
xt = 0.05 * (xl(2)-xl(1)) + xl(1);
yt = 0.90 * (yl(2)-yl(1)) + yl(1);
caption = sprintf('Outlet = %.3f * Inlet + %.3f', b(2), b(1));
text(xt, yt, caption, 'FontSize', 12, 'Color', 'k', 'FontWeight', 'bold');


%% Regulator Performance
setp = 400; % regulator set pressure, psi
SPE_DP = (press-setp)*SPE;
oxtank_spe = oxtank + SPE_DP;

figure
plot(time0(ind), oxtank_spe(ind));

%%
mdot_lox_2(1:ind(1)) = 0;
Q_ind = (ind(1)+50):(ind(end)-100); % taking only the relevant data
Q_lox = mdot_lox_2 /rho_lox; % volumetric flow of ox, ft3/s
Q_fuel = mdot_fuel_2 / rho_fuel; % volumetric flow of fuel, ft3/s
Q_tot = Q_lox + Q_fuel; % total flow, ft3/s
Q_scfm = Q_tot .* oxtank./ 14.7 .* 60;     % Total flow in standard cubic feet per minute, SCFM
%%
figure;
plot(time0(ind),Q_scfm(ind))
xlabel('Time (s)')
ylabel('Flow Rate SCFM')
ylim([20 100])
grid on
%%
Q_scfm_N2 = Q_scfm/2.56; % scaled for nitrogen flow rate for flow curve comparison
figure
pointsize=10;
scatter(Q_scfm_N2(Q_ind), oxtank_spe(Q_ind),pointsize,press(Q_ind),'x','LineWidth',5);
xlabel('Flow Rate (SCFM)'); ylabel('Outlet Pressure (psi)')
colormap jet
cb = colorbar();
ylim([200, 400])
xlim([10, 32])
grid on
h = colorbar;
ylabel(h, 'Inlet Pressure (psi)')

% Create scattered interpolant, to be used for future simulations
F_regoutlet = scatteredInterpolant(press(Q_ind),Q_scfm(Q_ind),oxtank_spe(Q_ind));

%% performance parameters
totimp = cumtrapz(time0(ind), Thrust(ind));
totimp = totimp(end);
totfuel = cumtrapz(time0(burnind), mdot_fuel_2(burnind));
totfuel= totfuel(end);
totlox = cumtrapz(time0(burnind), mdot_lox_2(burnind));
totlox = totlox(end);
totprop = totlox + totfuel;

mdot_lox_avg = mean(mdot_lox_2(burnind));
mdot_fuel_avg = mean(mdot_fuel_2(burnind));
mdot_avg = mean(mdot(burnind));
OF_avg = mean(OF(burnind));
CC_avg = mean(CC(burnind));
cstar_avg = mean(cstar(burnind));
cstar_eff_avg = mean(cstar_eff(burnind));
Thrust_avg = mean(Thrust(burnind));
Isp_avg = mean(Isp(burnind));

fprintf("Total Impulse: %.1f lbf-s \n", totimp)
fprintf("Ox Mass Flow Rate: %.3f lb/s \n", mdot_lox_avg)
fprintf("Fuel Mass Flow Rate: %.3f lb/s \n", mdot_fuel_avg)
fprintf("Mass Flow Rate: %.2f lb/s \n", mdot_avg)
fprintf("OF: %.2f \n", OF_avg)

fprintf("C-star: %.2f ft/s \n", cstar_avg)
fprintf("C-star Efficiency: %% %.2f \n", cstar_eff_avg)
fprintf("Isp: %.2f \n", Isp_avg)

fprintf("Chamber Pressure: %.1f psi \n", CC_avg)
fprintf("Thrust: %.1f lbf \n", Thrust_avg)

fprintf("Total Fuel: %.1f lbs \n", totfuel)
fprintf("Total Oxidizer Burned: %.1f lbs \n", totlox)
fprintf("Total Propellant Burned: %.1f lbs \n \n", totprop)


%% Chugging Analysis
avCC = smoothdata(CCRaw(burnind),'gaussian',10);
figure
plot(time0(burnind),CCRaw(burnind),time0(burnind),avCC)

%%
L = length(time0(burnind))*4; % increase length for better interpolation
tq = linspace(min(time0(burnind)), max(time0(burnind)), L); % create time array
inst = CCRaw(burnind)-avCC; % instbility, average CC subtracted out
CCq = interp1(time0(burnind),CCRaw(burnind),tq); % queried CC at sampling points
avCCq = interp1(time0(burnind),avCC,tq);  % queried average CC at sampling points
instq = CCq-avCCq;  % queried instability

% see how they match up
figure
plot(time0(burnind),inst,tq,instq);
legend('raw','interp')

T = mean(diff(tq)); % sampling period
Fs = 1/T; % Sampling Frequency
FF = fft(instq); % Fast Fourier Transform
P2 = abs(FF/L); % two-sided spectrum (includes negative frequencies)
P1 = 2*P2(1:L/2+1); % single-sided spectrum, x2 to account for all values
f = Fs*(0:(L/2))/L;  % create frequency array
figure
semilogx(f,smoothdata(P1,'gaussian',4),'LineWidth',2)
grid on
xlabel('Frequency, Hz')
ylabel('Amplitude')

