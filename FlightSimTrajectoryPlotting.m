%Flight Sim trajectory plotting
clear variables; close all; clc;
load('FlightSimApogeeSpread_T.mat');
rx_arr_T = rx_arr;
ry_arr_T = ry_arr;
tai_arr_T = tai_arr;
% load('FlightSimAOASpread_Tsub30.mat');
% aoa_arr_Tsub30 = aoa_arr;

rx_mean_T = mean(rx_arr_T(:,1:1000),1);
ry_mean_T = mean(ry_arr_T(:,1:1000),1);
aoa_std_T = std(rx_arr_T(:,1:1000),0,1);
tai_max_T = max(tai_arr_T);


aoa_mean_Tsub30 = mean(aoa_arr_Tsub30(:,1:1000),1);
aoa_std_Tsub30 = std(aoa_arr_Tsub30(:,1:1000),0,1);
t = t(1:1000);

sigma_count = 2;
curve_count = 2*sigma_count + 1;
curve_mean = sigma_count+1; % The curve number of the mean aoa curve

aoa_T = zeros(curve_count,length(t));
aoa_Tsub30 = zeros(curve_count,length(t));
for i = -sigma_count:sigma_count
    aoa_T(i+1+sigma_count,:) = rx_mean_T + i*aoa_std_T;
    aoa_Tsub30(i+1+sigma_count,:) = aoa_mean_Tsub30 + i*aoa_std_Tsub30;
end


%% Plotting
% Curves
curve_plots = gobjects(1,curve_count);
cmap = flipud(colormap(parula(2)));
hold on
for sigma_no = -sigma_count:sigma_count
    i = sigma_no+sigma_count+1;
    if sigma_no == 0 % Haven't gotten to increasing sigma
        curve_plots(i) = plot(t,aoa_T(i,:),'Color','red','LineWidth',2.5);
    elseif sigma_no < 0
        curve_plots(i) = patch([t,fliplr(t)],[aoa_T(i,:),fliplr(aoa_T(i+1,:))],cmap(abs(sigma_no),:),'linestyle','none','FaceAlpha',0.3);
    else
        curve_plots(i) = patch([t,fliplr(t)],[aoa_T(i-1,:),fliplr(aoa_T(i,:))],cmap(abs(sigma_no),:),'linestyle','none','FaceAlpha',0.3);
    end
end

%% Formatting
xlabel('Downrange E (mi)');
ylabel('Downrange S (mi)');

% Legend
for i = sigma_count:1
    curve_legend(i) = curve_plots(i + curve_mean);
    curve_string(i) = "\pm" + i + "\sigma";
end

curve_legend(end+1) = curve_plots(curve_mean);
curve_string(end+1) = "\mu";

legend(curve_legend,curve_string,'FontSize',12)