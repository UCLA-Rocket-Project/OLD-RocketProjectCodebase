function param = TargetVolumeParameter(map,loc,inorout,paramname)

height = size(map,1);

[h_loc, w_loc] = matsplit(loc);

if h_loc == 0 && w_loc == 0
    error('name could not be located in map')
elseif h_loc == 1 || h_loc == height
    error('name cannot be at bounds of map (since it should be a valve)')
elseif ~strcmp(inorout,'in') && ~strcmp(inorout,'out')
    error('inorout input must be ''in'' or ''out''')
end

if strcmp(inorout,'out')
    if ~isempty(map{h_loc+1,w_loc})
        param = map{h_loc+1,w_loc}.(sprintf('%s',paramname));
    else
        error('nothing found on outlet of name')
    end
elseif strcmp(inorout,'in')
    for dw = w_loc:-1:1
        if ~isempty(map{h_loc-1,dw})
            param = map{h_loc-1,dw}.(sprintf('%s',paramname));
            return
        end
    end
    error('nothing found on inlet of name')
end