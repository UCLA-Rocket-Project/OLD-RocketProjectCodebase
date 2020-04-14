function [T_arr, tb] = getThrust(T_xlsx, dt)
% Generates thrust array with uniform timesteps of dt using interpolation.
% T_data should be an xlsx file in the form 'FILENAME.xlsx'
% T_data should also start from t=0 and end at t=tb
% dt should be a power of 10

% Import data from file
T_data = readmatrix(T_xlsx);
t = T_data(:,1);
T = T_data(:,2);

% Get tb
tb = dt*floor(t(end)/dt); % no good way of rounding to nearest timestep (dt)

% Interpolate
tq = 0:dt:tb;
T_arr = interp1(t,T,tq);

end