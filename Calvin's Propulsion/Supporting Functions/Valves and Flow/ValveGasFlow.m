function q = ValveGasFlow(Cv,p_1,p_2,T_1,G_g,gam_g)
%   GasFlow: Calculate standard volumetric flow rate
%   Unit: std ft^3 (SCFM)
%   Cv = flow coefficient
%   p1 = upstream pressure in psia
%   p2 = downstream pressure in psia
%   G_g = gas specific gravity wrt air
%   T1 = upstream temperature in degR
%   gam_g = gas specific ratio

p_ratio = (2/(gam_g + 1))^(gam_g/(gam_g-1));

if (p_2 <= p_ratio*p_1)
    q = ValveChokedGasFlow(Cv,p_1,T_1,G_g);
else
    q = ValveUnchokedGasFlow(Cv,p_1,p_2,T_1,G_g);
end

    