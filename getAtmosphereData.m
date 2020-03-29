function [mu_X, sigma_all, no_time] = getAtmosphereData(Atm_xlsm)
% Reads FAR Atmosphere data from excel macro sheet. Non-loop function

atm_data = readmatrix(Atm_xlsm);

%% Get rid of NaN
% Empty spaces interpreted as NaN
[data_rows, data_cols] = size(atm_data);

i = 1; % used as index
j = 1; % used as counter
while j <= data_cols % Flip through manually
   col = atm_data(:,i); % Take single column
   is_nan = isnan(col); % Check for empty spaces
   if all(is_nan) % If empty column
       atm_data(:,i) = []; % Delete column
   else
       i = i + 1; % Only advance if column is not empty
   end
   j = j + 1;
end
[~, data_cols] = size(atm_data); % Update num columns 

%% Interpolate empty spots
for j = 1:data_cols
    for i = 1:data_rows
       if isnan(atm_data(i,j))
          no_empty = 1; % Count the empty cell
          i_sub = i;
          while isnan(atm_data(i_sub+1,j))
              no_empty = no_empty+1;
              i_sub = i_sub+1; % Final value is index for last empty cell
          end
          delta = (atm_data(i+no_empty,j)-atm_data(i-1,j)) / (no_empty+1);
          for k = i:i_sub
             atm_data(k,j) = (k - (i-1))*delta + atm_data(i-1,j);
          end
       end
    end    
end

%% For Each Time of Day: Get Data Array 
% (each column is a different day, each row is a different variable)
no_var = 3; % Temp, wind, P
no_day = data_cols/6; % Number of days that there is data for
no_time = 31; % 7 - 12 AM, 10 min increments

atm_var = zeros(no_var,no_day,no_time);
for time = 1:data_rows
    for data_col = 1:data_cols
        var = mod(data_col,6);
        day = floor((data_col-1)/6) + 1;
        switch var
            case 1
                atm_var(1,day,time) = atm_data(time,data_col); % (F) Temp
            case 4
                atm_var(2,day,time) = atm_data(time,data_col) * 1.46667; % (ft/s) wind
            case 0
                atm_var(3,day,time) = atm_data(time,data_col) * (14.696*144)/29.9212 ; % (lbf/ft^2) Pressure
        end
    end
end

%% Prep Covariance Matrices
% See Eq. 1 in https://en.wikipedia.org/wiki/Covariance_matrix

sigma_all = zeros(no_var,no_var,no_time); % Covariance matrices for all considered launch times
XX = zeros(no_var,no_var,no_day,no_time); % Unaveraged variable product
E_XX = zeros(no_var,no_var,no_time); % Averaged variable product
mu_X = zeros(no_var,no_time); % Averaged variables

for k = 1:no_time
    for j = 1:no_day
        XX(:,:,j,k) = atm_var(:,j,k)*transpose(atm_var(:,j,k)); % Calculate variable product for each set of data
    end
    E_XX(:,:,k) = mean(XX(:,:,:,k),3); % Mean of XX across all days for each time
    mu_X(:,k) = mean(atm_var(:,:,k),2); % Mean of each variable across all days for each time
    sigma_all(:,:,k) = E_XX(:,:,k) - mu_X(:,k)*transpose(mu_X(:,k));
end

