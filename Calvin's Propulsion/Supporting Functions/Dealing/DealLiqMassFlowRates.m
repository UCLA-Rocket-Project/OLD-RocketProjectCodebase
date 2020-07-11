function map = DealLiqMassFlowRates(map,loc,mdot)

%function param = TargetVolumeParameter(map,loc,inorout,paramname)

height = size(map,1);

[h_loc, w_loc] = matsplit(loc);

if h_loc == 0 && w_loc == 0
    error('name could not be located in map')
elseif h_loc == 1
    error('name cannot be at the top of the map (since it should be an injector)')
end

for dw = w_loc:-1:1
    if ~isempty(map{h_loc-1,dw})
        if (isa(map{h_loc-1,dw},'PropellantTank'))
            map{h_loc-1,dw}.vmdotliq = map{h_loc-1,dw}.vmdotliq + [0;0;mdot;0;0;-mdot];
        else
            error('invalid inlet object type')
        end
        return
    end
end
error('nothing found on inlet of name')