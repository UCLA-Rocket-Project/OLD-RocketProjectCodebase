function RecordStates(i,map)

height = size(map,1);
width = size(map,2);

for h = 1:1:height
    for w = 1:1:width
        if (~isempty(map{h,w}) && strcmp(map{h,w}.type,'volume'))
            map{h,w}.RecordStates(i);
        end
    end
end