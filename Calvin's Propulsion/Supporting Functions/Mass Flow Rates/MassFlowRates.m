function mdots = MassFlowRates(map,Xvals,A,dt,gas,comb_eff,combdata)

map = Resetmdots(map);

height = size(map,1);
width = size(map,2);

for h = 1:1:height
    for w = 1:1:width
        if (~isempty(map{h,w}) && ~strcmp(map{h,w}.type,'volume'))
            if strcmp(map{h,w}.phase,'gas')
                if isa(map{h,w},'Regulator')
                    loc_reg = [h,w];
                else
                    p_1 = TargetVolumeParameter(map,[h,w],'in','p');
                    T_1 = TargetVolumeParameter(map,[h,w],'in','T');
                    p_2 = TargetVolumeParameter(map,[h,w],'out','p');
                    mdot = map{h,w}.GasMassFlowRate(p_1,p_2,T_1,gas);
                    
                    map = DealGasMassFlowRates(map,[h,w],mdot);
                end
                    
            elseif strcmp(map{h,w}.phase,'liquid')
                if (strcmp(map{h,w}.type,'oxidizer'))
                    loc_ox = [h,w];
                elseif (strcmp(map{h,w}.type,'fuel'))
                    loc_fu = [h,w];
                end
            end
        end
    end
end

[mdot_liq_ox,mdot_liq_fu] = BipropInjectorMassFlowRate(map,loc_ox,loc_fu,comb_eff,combdata);
map = DealLiqMassFlowRates(map,loc_ox,mdot_liq_ox);
map = DealLiqMassFlowRates(map,loc_fu,mdot_liq_fu);

[h_reg,w_reg] = matsplit(loc_reg);
mdotreg = map{h_reg,w_reg}.GasMassFlowRate(loc_reg,map,Xvals,A,dt,gas);
map = DealGasMassFlowRates(map,loc_reg,mdotreg);

mdots = [];
volumes = VectorizeMapVolumes(map);
for j = 1:1:length(volumes)
    if ~isa(volumes(j),'CombustionChamber')
        volumes{j}.summdots;
        mdots = vertcat(mdots,volumes{j}.vmdot);
    end
end
