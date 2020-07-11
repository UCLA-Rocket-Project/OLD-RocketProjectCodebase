function map = MapInParallel(map,newobj,existingobj)
% MapInParallel: Add newobj in parallel with exisitingobj
% NOTE: Parallel objects shall never rejoin downstream

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
                if strcmp(map{h,w}.name,existingobj.name)
                    if w == width || isempty(map{h,w+1})
                        map{h,w+1} = newobj;
                        return
                    else
                        for dw = w+1:1:width
                            if isempty(map{h,w+1})
                                map{h,w+1+dw} = newobj;
                                return
                            end
                        end
                        map{h,dw+1} = newobj;
                        return
                    end
                end
            end
        end
        error('existing object input does not exist')
    end
end