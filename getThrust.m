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
len_data = length(t);
tb = dt*floor(t(len_data)/dt);
% len_arr = tb/dt + 1;

% Interpolate
tq = 0:dt:tb;
T_arr = interp1(t,T,tq); 


% Initialize T_arr
% T_arr = zeros(1,len_arr);
% T_arr(1) = T(1);

% Fix dt
% step_arr = 1;
% for i = 2:len_data
%     if (t(i)-t(i-1)) > dt
%         substeps = floor(t(i)/dt) - floor(t(i-1)/dt); % Find how many steps (for T_arr) needed for this i-value
%         m = T(
%         for j = 1:substeps
%             
%             step_arr = step_arr + 1; 
%         end
%     else
%         error('thrust data is too fine. add new method');
%     end
% end

end