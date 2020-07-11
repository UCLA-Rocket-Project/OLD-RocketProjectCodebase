function bool = OutOfPropellantListener(map)

height = size(map,1);
width = size(map,2);

for h = 1:1:height
    for w = 1:1:width
        if (~isempty(map{h,w}) && isa(map{h,w},'PropellantTank'))
            if map{h,w}.m_prop <= 0
                bool = true;
                sprintf('%s out of propellant',map{h,w}.name)
                return
            end
        end
    end
end
bool = false;