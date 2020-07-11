function map = MapInSeries(map,newobj,existingobj)
% MapInSeries: Add newobj in series with exisitingobj

height = size(map,1);
width = size(map,2);
switch nargin
    case 2
        if height ~= 0 && width ~= 0
            error('map size is nonzero. map must be new to add new object with no existing object as reference')
        else
            map{1,1} = newobj;
            return
        end
    case 3
        for h = 1:1:height
            for w = 1:1:width
                if ~isempty(map{h,w}) && strcmp(map{h,w}.name,existingobj.name)
                    if h == height || isempty(map{h+1,w})
                        map{h+1,w} = newobj;
                        return
                    else
                        error('next cell after exisiting object not empty')
                    end
                end
            end
        end
        error('existing object input does not exist')
    end
end