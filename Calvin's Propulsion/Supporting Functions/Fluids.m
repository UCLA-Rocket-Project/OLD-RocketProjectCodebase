function fluids = Fluids()

He = Gas;
He.R = 12421; %[ft lbf/slug °R]
He.gam = 1.66;
He.G = 0.138;
He.c_v = He.R/(He.gam - 1);
He.c_p = He.R*He.gam/(He.gam - 1);

% He = Gas;
% He.R = 1775; %[ft lbf/slug °R]
% He.gam = 1.4;
% He.G = 1;
% He.c_v = He.R/(He.gam - 1);
% He.c_p = He.R*He.gam/(He.gam - 1);

LOx = Liquid;
LOx.rho = 71.17/32.2;

Eth = Liquid;
Eth.rho = 50.11/32.2;

fluids = struct('He',He,'LOx',LOx,'Eth',Eth);