%% load data
data = csvread('Static Fire 2-23-20.csv', 1, 0);
% data = xlsread('hotfire2.xlsx');
 timeRaw = data(:, 1);
% oxraw = data(:, 6);
% ox = 244.95*oxraw - 239.92;
% thrust = data(:, 9);
% press = data(:, 11);
% anothermanifold = data(:, 10);
% cc = data(:, 12);
% manifold = data(:, 13);
oxman = data(:, 3);
fuelman = data(:, 5);
cc = data(:, 7);
ox = data(:, 9);
press = data(:, 11);
thrust = data(:, 13);
thrusto = thrust - 22;
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

%% Calibrate data
time0 = 1905448;
time = (timeRaw-time0)/1000;
ind = 4000:9000;

load 'Trial4'

%%
figure(1)
hold on
plot(time(ind), thrusto(ind), 'LineWidth', 3)
plot(t(1:i),map{7,1}.data.FT(1:i),'linewidth',3)
hold off
title('2/20/20 Thrust minus LC offset')
grid on
grid minor
xlim([-1 18]);
%%
figure(2)
hold on
plot(time(ind), ox(ind), 'LineWidth', 4);
plot(t(1:i),map{5,1}.data.p(1:i),'linewidth',lw)
hold off
title('Oxygen Tank')
xlabel('Time [s]')
ylabel('Tank Pressure [psi]')
grid on
grid minor
xlim([-10, 18]);
%%
At = 1.77; % in2
C_tau = 1.4;
CCtheory = thrusto / (At*C_tau);
thrusttheory = cc*At*C_tau;
figure(4)
plot(time(ind), thrust(ind), time(ind), thrusto(ind), time(ind), thrusttheory(ind), 'LineWidth', 3)
legend('Load Cell', 'Load Cell minus 29.6 lbs', 'Calculated from chamber pressure')
xlim([-1, 18])
ylim([0, 700]);
set(gca, 'FontSize', 17, 'FontWeight', 'bold')
grid on
grid minor
ylabel('Pressure [psi]')
xlabel('Time [s]')
title('Thrust')
figure(3)
plot(time(ind), cc(ind), time(ind), CCtheory(ind), 'LineWidth', 3)
legend('PT data', 'Calculated from Thrust minus 29.6 pounds')
xlim([-1, 20])
grid on
grid minor
ylabel('Pressure [psi]')
xlabel('Time [s]')
title('Chamber Pressure')