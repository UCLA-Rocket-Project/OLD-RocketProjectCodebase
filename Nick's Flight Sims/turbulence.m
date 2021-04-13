function [vwgturb,vwg] = turbulence(vwg,turb,z,s,nufs,vwg20,vwgturbc,dt,it)
    vwgturbc1 = vwgturbc(1,1);
    vwgturbc2 = vwgturbc(1,2);
    if z>=21 && z<=1000
        Lw = z;
        Lu = z/(0.177+0.000823*z)^1.2;
        sigw = 0.1*vwg20;
        sigu = sigw/(0.177+0.000823*z)^0.4;
    elseif z>=2000
        Lw = 1750;
        Lu = 1750;
        if turb == 1  % Change to the high altitude chart later
            sigw = 5;
        elseif turb == 2
            sigw = 10;
        elseif turb == 3
            sigw = 15;
        end
        sigu = sigw;
    elseif z>1000 && z<2000
        Lw = 1000+(z-1000)*750/1000;
        Lu = Lw;
        if turb == 1  % Change to the high altitude chart later
            sigw = 0.1*vwg20+(z-1000)*(5-0.1*vwg20)/1000;
        elseif turb == 2
            sigw = 0.1*vwg20+(z-1000)*(10-0.1*vwg20)/1000;
        elseif turb == 3
            sigw = 0.1*vwg20+(z-1000)*(15-0.1*vwg20)/1000;
        end
        sigu = sigw;
    end
    %         Gun(i,:) = sigu(i)*sqrt(2*Lu(i)/pi/s(i))*[0 1];
    %         Gud(i,:) = [Lu(i)/s(i) 1];
    %         [AuT(i),BuT(i),CuT(i),DuT(i)] = tf2ss(Gun(i,:),Gud(i,:));
    %         vwgturb(i,1) = CuT(i)*xTurb(i)+DuT(i)*nufs(i,1);
    %         xTurb(i+1) = AuT(i)*xTurb(i)+BuT(i)*nufs(i,1);
    % %         vwgturb(i,1) = mean(filter(Gun(i,:),Gud(i,:),nufs(:,i)));
    %         Gwn(i,:) = sigw(i)*sqrt(2*Lw(i)/pi/s(i))*[2*sqrt(3)*Lw(i)/s(i) 1];
    %         Gwd(i,:) = [4*(Lw(i)^2)/(s(i)^2) 4*Lw(i)/(s(i)) 1];
    %         [AwT,BwT,CwT,DwT] = tf2ss(Gwn(i,:),Gwd(i,:));
    %         vwgturb(i,2) = CwT*zTurb(i,:).'+DwT*nufs(i,2);
    %         zTurb(i+1,:) = AwT*zTurb(i,:).'+BwT*nufs(i,2);
    % %         vwgturb(i,2) = mean(filter(Gwn(i,:),Gwd(i,:),nufs(:,i)));
    vwgturb(1,1) = (1-s*(dt/100)/Lu)*vwgturbc1+...
        sqrt(2*s*(dt/100)/Lu)*sigu/std(nufs(1,:))*nufs(1,it);
    vwgturb(1,2) = (1-s*(dt/100)/Lu)*vwgturbc2+...
        sqrt(2*s*(dt/100)/Lu)*sigu/std(nufs(2,:))*nufs(2,it);
    vwg(1:2) = vwg(1:2) + vwgturb;
end