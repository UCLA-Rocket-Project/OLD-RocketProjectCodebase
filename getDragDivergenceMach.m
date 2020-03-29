function [Ma_dd] = getDragDivergenceMach(RAS)
% The definition of drag divergence mach is Ma at which
% dCD/dMa = 0.1 (i.e. Mach where drag increases abruptly)
Ma = transpose(RAS(:,1));
CD = transpose(RAS(:,3));

Ma_dd = 0;
dMa = 0.01;
for i = 2:length(Ma)-1 % Assume no wackiness for small i
    if (CD(i+1)-CD(i-1))/(2*dMa) > 0.1 && Ma(i)>0.4 % Assume increments of Ma = 0.01, linearizing
        
        dCD2 = (CD(i+1)-CD(i-1))/(2*dMa); %(CD(i+1)-2*CD(i)+CD(i-1))/dMa^2;
        dCD1 = (CD(i)-CD(i-2))/(2*dMa); %(CD(i+1)-2*CD(i)+CD(i-1))/dMa^2;
        m = (dCD2-dCD1)/dMa; % linearize second derivative
        b = dCD2 - m*Ma(i);
        assert(abs(m*Ma(i-1)+b - dCD1) < 10^-3);
        
        Ma_dd = (0.1-b)/m; % find Ma_dd linearly
        
        break
    end
end