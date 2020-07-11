function DefineGoverningEquations(map,gas)

height = size(map,1);
width = size(map,2);

for h = 1:1:height
    for w = 1:1:width
        if (~isempty(map{h,w}) && (isa(map{h,w},'PropellantTank') || isa(map{h,w},'Tubing')))
            for dw = w:-1:1
                if ~isempty(map{h-2,dw})
                    map{h,w}.DefineGoverningEquations(map{h-2,dw},gas);
                    break
                end
            end
        elseif isa(map{h,w},'PressurantTank')
            map{h,w}.DefineGoverningEquations(gas);
        end
    end
end
