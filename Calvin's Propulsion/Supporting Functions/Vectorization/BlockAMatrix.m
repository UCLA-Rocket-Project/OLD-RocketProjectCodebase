function A = BlockAMatrix(volumes)
% DefineXvals: Defines and outputs numerical vertical vector of object state values
A = [];
for i = 1:1:length(volumes)
    A = blkdiag(A,volumes{i}.A);
end

if isempty(A)
    error('A did not update; nothing')
end