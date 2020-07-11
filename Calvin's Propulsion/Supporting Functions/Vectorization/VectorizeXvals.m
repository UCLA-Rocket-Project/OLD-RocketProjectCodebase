function Xvals = VectorizeXvals(volumes)
% DefineXvals: Defines and outputs numerical vertical vector of object state values
Xvals = [];
for i = 1:1:length(volumes)
    volumes{i}.DefineXvals();
    Xvals = horzcat(Xvals,volumes{i}.Xvals);
end

Xvals = transpose(Xvals);