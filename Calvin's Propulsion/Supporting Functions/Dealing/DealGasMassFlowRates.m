function map = DealGasMassFlowRates(map,loc,mdot)

%function param = TargetVolumeParameter(map,loc,inorout,paramname)

height = size(map,1);

[h_loc, w_loc] = matsplit(loc);

if h_loc == 0 && w_loc == 0
    error('name could not be located in map')
elseif h_loc == 1 || h_loc == height
    error('name cannot be at bounds of map (since it should be a valve)')
end

if ~isempty(map{h_loc+1,w_loc})
    %param = map{h_loc+1,w_loc}.(sprintf('%s',paramname));
    if (isa(map{h_loc+1,w_loc},'Tubing'))
        map{h_loc+1,w_loc}.vmdotin = map{h_loc+1,w_loc}.vmdotin + [0;mdot;0;0];
        map{h_loc+1,w_loc}.vmdotspec = map{h_loc+1,w_loc}.vmdotspec + map{h_loc+1,w_loc}.specfactor.*[0;0;0;mdot];
    elseif (isa(map{h_loc+1,w_loc},'PropellantTank'))
        map{h_loc+1,w_loc}.vmdotin = map{h_loc+1,w_loc}.vmdotin + [0;mdot;0;0;mdot;0];
    else
        error('invalid outlet object type')
    end
else
    error('nothing found on outlet of name')
end

for dw = w_loc:-1:1
    if ~isempty(map{h_loc-1,dw})
        if (isa(map{h_loc-1,dw},'Tubing'))
            map{h_loc-1,dw}.vmdotout = map{h_loc-1,dw}.vmdotout + [0;-mdot;0;-mdot];
        elseif (isa(map{h_loc-1,dw},'PressurantTank'))
            map{h_loc-1,dw}.vmdotout = map{h_loc-1,dw}.vmdotout + [0;-mdot;0;-mdot];
        elseif dw == 1
            error('invalid inlet object type')
        end
        return
    end
end
error('nothing found on inlet of name')