function DealXvals(volumes,Xvals)
k = 1;
for i = 1:1:length(volumes)
    if isa(volumes{i},'PressurantTank')
        volumes{i}.p = Xvals(k);
        volumes{i}.T = Xvals(k+1);
        volumes{i}.rho = Xvals(k+2);
        volumes{i}.m = Xvals(k+3);
        k = k + 4;
    elseif isa(volumes{i},'Tubing')
        volumes{i}.p = Xvals(k);
        volumes{i}.T = Xvals(k+1);
        volumes{i}.rho = Xvals(k+2);
        volumes{i}.m = Xvals(k+3);
        k = k + 4;
    elseif isa(volumes{i},'PropellantTank')
        volumes{i}.p = Xvals(k);
        volumes{i}.T = Xvals(k+1);
        volumes{i}.rho = Xvals(k+2);
        volumes{i}.m = Xvals(k+3);
        volumes{i}.V = Xvals(k+4);
        volumes{i}.m_prop = Xvals(k+5);
        k = k + 6;
    else
        error('invalid object type')
    end
end