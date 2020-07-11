function map = Resetmdots(map)

height = size(map,1);
width = size(map,2);

for h = 1:1:height
    for w = 1:1:width
        if (~isempty(map{h,w}) && (isa(map{h,w},'PropellantTank') || isa(map{h,w},'PressurantTank') || isa(map{h,w},'Tubing')))
            map{h,w}.Resetmdots;
        end
    end
end