function objv = VectorizeMapVolumes(map)

height = size(map,1);
width = size(map,2);

objv = {};

k = 1;
for h = 1:1:height
    for w = 1:1:width
        if ~isempty(map{h,w}) && strcmp(map{h,w}.type,'volume') && ~isa(map{h,w},'CombustionChamber')
            objv{k} = map{h,w};
            k = k + 1;
        end
    end
end
