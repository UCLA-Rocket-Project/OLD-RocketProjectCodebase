function [mdot_liq_ox,mdot_liq_fu] = BipropInjectorMassFlowRate(map,loc_ox,loc_fu,comb_eff,comb_data)

[h_ox,w_ox] = matsplit(loc_ox);
[h_fu,w_fu] = matsplit(loc_fu);

cstar_ideal = lininterp1(comb_data.OF_ref,comb_data.cstar_ref,map{h_fu+1,w_fu}.OF);          % ideal characteristic velocity [ft/s]
map{h_fu+1,w_fu}.cstar = cstar_ideal*comb_eff*.01;                       % real characteristic velocity [ft/s]

n = 25;

vp_c_fu = linspace(0,map{h_fu-1,w_fu}.p,n);
vp_c_ox = linspace(0,map{h_ox-1,w_ox}.p,n);

vdP_fu = map{h_fu-1,w_fu}.p - vp_c_fu;
vdP_ox = map{h_ox-1,w_ox}.p - vp_c_ox;

vmdot_fu = map{h_fu,w_fu}.CdA*sqrt(2*map{h_fu-1,w_fu}.rho_prop*vdP_fu/144);
vmdot_ox = map{h_ox,w_ox}.CdA*sqrt(2*map{h_ox-1,w_ox}.rho_prop*vdP_ox/144);
vmdot_tot = vmdot_fu + vmdot_ox;

vp_c_noz = vmdot_tot*map{h_fu+1,w_fu}.cstar/map{h_fu+1,w_fu}.A_t;
map{h_fu+1,w_fu}.p = lininterp1(vp_c_noz - vp_c_fu,vp_c_fu,0);

mdot_liq_fu = lininterp1(vp_c_noz-vp_c_fu,vmdot_fu,0);
mdot_liq_ox = lininterp1(vp_c_noz-vp_c_fu,vmdot_ox,0);
map{h_fu+1,w_fu}.mdot = mdot_liq_fu + mdot_liq_ox;

map{h_fu+1,w_fu}.OF = mdot_liq_ox/mdot_liq_fu;

map{h_fu+1,w_fu}.Thrust;