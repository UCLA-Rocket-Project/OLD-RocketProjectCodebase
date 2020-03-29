clear variables; clc; close all;

z0 = 0; zf = 50000; dz = 0.1;
steps = ((zf-z0)/dz)+1;
z = linspace(z0,zf,steps);

rail = 60; %60 ft rail

wind = zeros(4,steps);
for j = 1:3
    for i = 1:steps
        if z(i) <= 2000*3.28084
            wind(j,i) = windProfileLog((j*5)*1.46667,rail,z(i));
        elseif z(i) > 2000*3.28084
            wind(j,i) = windProfilePower((j*5)*1.46667,rail,z(i));
        end
    end
end
    
hold on;
plot1 = plot(z,wind(1,:)/1.46667,'LineWidth',3);
plot2 = plot(z,wind(2,:)/1.46667,'LineWidth',3); 
plot3 = plot(z,wind(3,:)/1.46667,'LineWidth',3);

xlabel('Altitude (ft)');
xlim([0 50000]);
ylabel('Wind Velocity (mph)');
title({'Flight Simulation Wind Model'})
lgd = legend([plot3 plot2 plot1],{'15 mph','10 mph','5 mph'});
legend('Location','northwest');
title(lgd,'Off-Rail Wind Velocity');
grid on