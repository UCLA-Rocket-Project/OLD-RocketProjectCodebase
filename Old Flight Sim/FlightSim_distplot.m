%Flight Sim trajectory plotting
clear variables; close all; clc;
load('FlightSim_dist.mat');
img_far = imread('FAR h12xw20.png');
img_h_mi = 12; img_w_mi = 20;
img_h = size(img_far,1); img_w = size(img_far,2);

radius_mean = mean(dist_arr);
radius_std = std(dist_arr);
sigma_count = 3;
circle_count = sigma_count + 1;
circle_count_n = 0;
for i = 1:sigma_count
    if (radius_mean - i*radius_std) >= 0
        circle_count = circle_count + 1;
        circle_count_n = circle_count_n + 1;
    else
        break
    end
end
circle_mean = circle_count_n+1; % The circle number of the mean landing circle

center = zeros(circle_count,2);
radius = zeros(circle_count,1);
for i = -circle_count_n:sigma_count
    radius(i+circle_count_n+1) = radius_mean + i*radius_std;
end

%center = center + [img_w/2 - 2, img_h/2 + 0]; %(-2,0) for 12x20; 
%radius = radius * img_h/img_h_mi;


%% Plotting
% Google Maps image
figure;
hold on;
% RI = imref2d(size(img_far));
% RI.XWorldLimits = [-img_w_mi/2 img_w_mi/2];
% RI.YWorldLimits = [-img_h_mi/2 img_h_mi/2];
% imshow(img_far,RI);

% Circles
t = linspace(0,2*pi);
circle_plots = gobjects(1,circle_count);
cmap = flipud(colormap(parula(3)));
for i = 1:circle_count % These are outer/enclosing circles
    if i == 1
        sigma_no = circle_count_n+1; % Innermost circle
        x0 = center(i,1); y0 = center(i,2);    % circles center
        r = radius(i);
        circle_plots(i) = patch([x0+r*cos(t)],[y0+r*sin(t)],cmap(sigma_no,:),'linestyle','-','FaceAlpha',0.3,'EdgeColor','black');
    else
        if i <= circle_mean % Haven't gotten to increasing sigma
            sigma_no = abs(i-circle_mean-1);
        else
            sigma_no = i - circle_mean;
        end
        x0 = center(i,1); y0 = center(i,2);    % circles center
        rin = radius(i-1); rout = radius(i);   % radii sizes
        circle_plots(i) = patch([x0+rout*cos(t),x0+rin*cos(t)],[y0+rout*sin(t),y0+rin*sin(t)],cmap(sigma_no,:),'linestyle','-','LineWidth',1.5,'FaceAlpha',0.3,'EdgeColor','black');
    end
end
circle_plots(i+1) = viscircles(center(circle_mean,:),radius(circle_mean),'Color','red','LineWidth',2.5);

%% Formatting
xlabel('Downrange E (mi)');
ylabel('Downrange S (mi)');

% Legend
for i = 1:sigma_count
    circle_legend(i) = circle_plots(circle_count_n + i + 1);
    circle_string(i) = "\pm" + i + "\sigma";
end
circle_legend = flipud(transpose(circle_legend));
circle_string = flipud(transpose(circle_string));

circle_legend(end+1) = circle_plots(end);
circle_string(end+1) = "\mu = " + round(radius_mean,2) + "mi";

legend(circle_legend,circle_string,'FontSize',12)