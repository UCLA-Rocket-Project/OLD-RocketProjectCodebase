function [u] = windProfilePower(groundWind,rail,z)
%Gets wind speed for z(altitude)>2000m.
%https://en.wikipedia.org/wiki/Wind_profile_power_law
%z is in ft, U is array in ft/s

z_ref = 2000*3.28084; %convert to ft
u_ref = windProfileLog(groundWind,rail,z_ref); %in ft/s
a = 1/7; %exponent in power law for neutral stability conditions

u = u_ref*(z/z_ref)^a; %ft/s