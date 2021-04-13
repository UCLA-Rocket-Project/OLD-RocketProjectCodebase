function [vwg,dpdd] = windspeed(vwgd,rho,P2,P1,H,z,h,vz,ustar,k,z1,LMBL,hv,fv,dpddp)
    if H == h
        vwg = vwgd;
        dpdd = 0;
    elseif z<=hv % Implement geostrophic patterns at z > hv
%         if atm_stability == 1
            vwg = ustar/k*(log(z/z1)+z/LMBL-z/hv*(z/2/LMBL));
%         elseif atm_stability == 2
%             vwg(i,1) = ustar(i)/k*(log(z/z1)+4.7*z/Lv*(1-z/2/hv)+z/LMBL-z/hv*(z/2/LMBL));
%         elseif atm_stability == 3
%             vx = (1-12*z/Lv)^(1/3);
%             fzl = 1.5*log((1+vx+vx^2)/3)-sqrt(3)*atan((1+2*vx)/sqrt(3))+pi/sqrt(3);
%             vwg(i,1) = ustar(i)/k*(log(z/z1)-fzl+z/LMBL-z/hv*(z/2/LMBL));
%         end
        if vz>=0
            dpdd = vwg*rho*fv;
%             dpdxi = i;
        end
    elseif z>hv
        if vwgd == 0
            dpdd = 0;
        else
            dpdd = ((P2/P1+1)/2)*dpddp;
%             dpdx(i) = dpdx(i-1);
%             dpdx(i) = d(P(i)/P(i-1))*dpdx(i-1);
        end
        vwg = dpdd/fv/rho;
    end
end