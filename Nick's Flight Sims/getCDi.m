function [CDab,CDaf] = getCDi(aoa,aoadd,aoadn,Rl,Ln,Afp,dfa,kfb,kbf,Afe,dist,RL,d)
    if abs(aoa)<20*pi/180 && abs(aoa)>=4*pi/180
        aoa1 = 2*floor((abs(aoa)*180/pi)/2);
        aoa2 = aoa1+2;
        i1 = aoa1/2-1;
        i2 = i1+1;
        del = aoadd(i2)-(aoadd(i2)-aoadd(i1))*(aoa2-abs(aoa)*180/pi)/(aoa2-aoa1);
        ne = aoadn(i2)-(aoadn(i2)-aoadn(i1))*(aoa2-abs(aoa)*180/pi)/(aoa2-aoa1);
    elseif abs(aoa)<=4*pi/180
%         del = aoadd(1)+(aoadd(2)-aoadd(1))*(abs(aoa(i))*180/pi)/2;
%         ne = aoadn(1)+(aoadn(2)-aoadn(1))*(abs(aoa(i))*180/pi)/2;
        del = aoadd(1);
        ne = aoadn(1);
    else
        del = aoadd(9);
        ne = aoadn(9);
    end
    
    if abs(aoa)<=20*pi/180 %20*pi/180
        CDab = 2*del*(aoa)^2+3.6*ne*(1.36*Rl-0.55*Ln)*abs(aoa)^3/(pi*(d/12));
        CDaf = (aoa)^2*(1.2*Afp*4/(pi*dfa^2)+3.12*(kfb+kbf-1)*Afe*4/(pi*dfa^2));
    else
        CDab = 2*del*(20*pi/180)^2+3.6*ne*(1.36*Rl-0.55*Ln)*(20*pi/180)^3/(pi*(d/12));
        CDaf = (20*pi/180)^2*(1.2*Afp*4/(pi*dfa^2)+3.12*(kfb+kbf-1)*Afe*4/(pi*dfa^2));
    end
    if dist < RL
        CDab = 0;
        CDaf = 0;
    end
%     CDab = 0;
%     CDaf = 0;
end