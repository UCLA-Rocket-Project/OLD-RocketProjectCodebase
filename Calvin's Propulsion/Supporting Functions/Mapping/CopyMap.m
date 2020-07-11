function newmap = CopyMap(map)
height = size(map,1);
width = size(map,2);

newmap = cell(height,width);

for h = 1:1:height
    for w = 1:1:width
        if (~isempty(map{h,w}))
            newmap{h,w} = copy(map{h,w});
        end
    end
end