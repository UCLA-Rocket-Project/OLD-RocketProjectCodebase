% Problem 1 Mod Solution:
% One-Dimensional Flight Modeling
clc; clear variables; close all;

%% Initial values
% Use values from Endurance
% Add initial/burn/final times, timestep, total timesteps
% Add wet/dry weight, gravity, average thrust
% Add initial conditions (a, v, z)

t_i = 0; t_b = 16; t_f = 160; dt = 0.01;
t = [t_i:dt:t_f];
t_steps = length(t);

W_wet = -139.80; % weight points down
W_dry = -91.53;
W_dot = (W_dry - W_wet)/t_b;
g = 32.174;

T_avg = 580;

v_i = 0; z_i = 0;

%% Initialize arrays
% We want to track:
% time
% weight, thrust, total force
% acceleration, velocity, position

W = zeros(1,t_steps);
m = zeros(1,t_steps);
T = zeros(1,t_steps);
F = zeros(1,t_steps);

atm_rho = zeros(1,t_steps);
atm_P = zeros(1,t_steps);
atm_T = zeros(1,t_steps);

a = zeros(1,t_steps);
v = zeros(1,t_steps);
z = zeros(1,t_steps);

%% Fill Initial Conditions
    W(1) = W_wet;
    m(1) = abs(W(1)) / g;
    T(1) = T_avg;
    
    v(1) = v_i;
    z(1) = z_i;
    
%% Simulate
% Follow pseudocode from slides

for i = 1:t_steps-1
    
    %% Current State Calculations  
    % Newton's 2nd Law
    F(i) = T(i) + W(i);
    a(i) = F(i) / m(i);
    
    %% Calculate Future State
        if t(i+1) <= t_b
            W(i+1) = W(i) + W_dot*dt;
            m(i+1) = abs(W(i+1)) / g;
            T(i+1) = T_avg;
        else
            W(i+1) = W_dry;
            m(i+1) = abs(W(i+1)) / g;
            T(i+1) = 0;
        end

        % Euler method
        v(i+1) = v(i) + a(i)*dt;
        z(i+1) = z(i) + v(i)*dt;
    
    %% End Condition
    if z(i) <= 0 && v(i) < 0
        break
    end
end

%% Results
% Find apogee, max velocity, max acceleration
% as well as the times at which they happen

[z_max, i_z_max] = max(z);
[v_max, i_v_max] = max(v);
[a_max, i_a_max] = max(a);

t_z_max = t(i_z_max);
t_v_max = t(i_v_max);
t_a_max = t(i_a_max);

%% Plotting
% Plot altitude vs. time
figure(1)
plot(t,z)
xlabel('Time (s)')
ylabel('Altitude (ft)')

% Plot velocity vs. time
figure(2)
plot(t,v)
xlabel('Time (s)')
ylabel('Vertical Velocity (ft/s)')

% Plot acceleration vs. time
figure(3)
plot(t,a)
xlabel('Time (s)')
ylabel('Vertical Acceleration (ft/s^{2})')