function [u] = windProfileLog(groundWind,rail,z)

%Convert to meters
groundWind = groundWind/3.28084;
rail = rail/3.28084;
z = z/3.28084;

%Initial units in metric
k = 0.41; %Von Karman constant
d = 0; %zero plane displacement, estimating for FAR site
z0 = 0.005; %roughness length for desert (usually 0.001 - 0.005 m)

u_star = groundWind*k/log((rail-d)/z0); %friction velocity, approx. for Mojave
% Assume weather data taken at rail length above ground

u = (u_star/k)*log((z-d)/z0); %log wind profile
u = u*3.28084; %convert to ft/s