%% Clear Contamination From Past Runs

%clc; 
clear variables; 
% clear variables Frfs 
% clearvars -except zzz
close all;
format short

%% Initial values & Initialize arrays

load('InputParams_6DOF');

%% Simulate

i = 1;
for t = ti:dt:tf
    %% There must be forces on rocket at all times
    if i>length(F)
%         disp(i)
        break
    end
    
    if t>tb
%         cdv = mach_cd(:,7);
    end
    
    %% Atmosphere Model
    
    z = x(i,3)-h;
    zs(i) = z;
    H = x(i,3);
    Hs(i) = H;
    
    % End Condition
    
    if z < 0 && t > tb && v(i,3) < 0
        tff = t;
        itf = i;
%     disp(i+2)
        break
    end
    
    % % % %
    
    if i==1
        tempr0 = 59 - 0.00356*H;
        H0 = H;
        Pr0 = 2116*((Temp0+459.67)/518.6)^5.256;
    end
    [temp(i),P(i),rho(i)] = atmosphere(H,Temp0,P0,tempr0,Pr0,H0,kTemp,R);
    
    hr(i) = hr0;
    rhod(i) = rho(i);
    rho(i) = rhod(i)*(1+hr(i))/(1+1.609*hr(i));
    
    %% Wind Speed Model

    if i~=1
        dpdxp = dpdx(i-1);
        dpdyp = dpdy(i-1);
        P1 = P(i-1);
    else
        dpdxp = 0;
        dpdyp = 0;
        P1 = P(i);
    end
    [vwg(i,1),dpdx(i)] = windspeed(vwgx,rho(i),P(i),P1,H,z,h,v(i,3),...
    ustar(i),k,z1,LMBL,hv,fv,dpdxp);
%     vwg(i,1) = -vwg(i,1);
    [vwg(i,2),dpdy(i)] = windspeed(vwgy,rho(i),P(i),P1,H,z,h,v(i,3),...
    vstar(i),k,z1,LMBL2,hv2,fv,dpdyp);
%     vwg(i,2) = -vwg(i,2);

    %% Gravity Model
    
    %g = (glat-3.086*(10^(-6))*(H/3.28))*3.28;
    g = GM/(re+H)^2;
    m = w(i)/g;
    
    %% Wind Turbulence Model
    
    if z >= 19 && z <= 21 % or 22 or 23
        vwg20 = vwg(i,1);
    end
    
    if turb ~= 0 && z > 21
        [vwgturb(i+1,:),vwg(i,:)] = turbulence(vwg(i,:),turb,z,s(i),nufs,vwg20,vwgturb(i,:),dt,i);
    end
    
    %% Wind Gusting Model
    
    vgust(i) = windgust(gust,t,tgust,dtgust,dtrise);
    vwg(i,1) = vwg(i,1) + vgust(i);
    
    %% Find speed relative to wind and angle of attack of the center of mass
    
    vrw(i,:) = v(i,:)-vwg(i,:);
    aoacm(i) = acos(dot((vrw(i,:)/norm(vrw(i,:))),rA(:,i)));
    if v(i,3)<0
        aoacm(i) = 0;
    end
    
    %% Find Mach Number and Set Compressibility Correction Factors
    
    vs(i) = (gamma*R*(temp(i)+459.69))^(1/2);
%     mach(i) = round(norm(v(i,:))/vs(i),2);
    mach(i) = norm(vrw(i,:))/vs(i);
    
    if mach(i)>0.9 && mach(i)<1.1
        beta = sqrt(1-0.9^2);
        CNaf = CNa(1,4)/beta;
%         CNaf = 2*pi*Afin*((Afe/2)/Aref)/(2+sqrt(4+beta*Afin/cos(thetaLE)));
    elseif mach(i) <= 0.9
        beta = sqrt(1-mach(i)^2);
        CNaf = CNa(1,4)/beta;
%         CNaf = 2*pi*Afin*((Afe/2)/Aref)/(2+sqrt(4+beta*Afin/cos(thetaLE)));
    else
        beta = sqrt(mach(i)^2-1);
        CNaf = CNa(1,4)/beta;
%         C3 = (gamma*mach(i)^4+(mach(i)^2-2)^2)/(2*(mach(i)^2-1)^1.5);
%         Aprime = 2*ft/3/vs(i);
%         CNaf = 4/beta*(1-(1-C3*Aprime)/2/beta/Afin);
%         CNaf = K_wb*(Sfin/Aref)*CNaf;
        nosebodyCNa = msCNa*mach(i) + bsCNa;
        CNan = 0.2*mach(i) + 2.03;
        CNab = nosebodyCNa - CNan;
    end
    betas(i) = beta;
    
    if mach(i)>1.1 && mach(i)<1.15
%         beta = sqrt(mach(i)^2-1);
    end
    
    %% Calculate Center of Pressure
    
%     CNaf = CNa(1,4)*beta;
    if abs(aoacm(i))<20*pi/180 && mach(i)<1.2
        CNab = 4*K/pi*Ap/(d/12)^2*abs(aoacm(i));
    elseif mach(i)<1.2
        CNab = 4*K/pi*Ap/(d/12)^2*20*pi/180;
    end
    %CNa = CNan+CNat+CNaf+CNab;
%     if mach(i)<1.2
        CNa(i,:) = [CNa(1,1)/beta,CNab/beta,CNa(1,3)/beta,CNaf];
%     else
%         CNa(i,:) = [CNan,CNab,CNa(1,3)/beta,CNaf];
%     end
%     CNa(i,:) = [CNa(1,1),CNab,CNa(1,3),CNa(1,4)];
%     CNa(i,:) = [CNa(1,1)/beta,CNab/beta,CNa(1,3)/beta,CNaf];

    %CP8(i) = sum(CNa(i,:).*CPi(1,:))/sum(CNa(i,:));
%     CPi(i,:) = CPi(1,:)/beta;
    if mach(i)>0.5
        CPi(i,4) = CPi(1,4)+(CPiFM2-CPiFM0)/(2-0.5)*(mach(i)-0.5);
    else
        CPi(i,4) = CPi(1,4);
    end
        
    CPi(i,1:3) = CPi(1,1:3);
%     CPi(i,:) = CPi(1,:);
%     CP8(i) = sum(abs(CNa(i,:)).*CPi(i,:))/sum(abs(CNa(i,:))); % Why abs?
    CP8(i) = sum(CNa(i,:).*CPi(i,:))/sum(CNa(i,:)); % Why abs?
%     CP(i) = CP8(i);
%     CPb(i) = sum(CNa(i,1:3).*CPi(i,1:3))/sum(CNa(i,1:3));
%     CP8(i) = CP8(i)+sin(abs(aoa(i)))*(CPcut-CP(i));

    %CP8(i) = CP8(i)/beta;
    if CP8(i)>Rl
        CP8(i)=Rl;
    elseif CP8(i)<0
        CP8(i)=0;
    end
    
    %% Thrust Model
    
    if t>tb
        Tt(i) = 0;
    end
    if i==2
        Tt(i) = 2*max(Tt);
    end
    if t<=tb 
%         tA(1:3,i) = [cos(tma(2)),0,sin(tma(2));0,1,0;-sin(tma(2)),0,cos(tma(2))]*...
%             [1,0,0;0,cos(tma(1)),-sin(tma(1));0,sin(tma(1)),cos(tma(1))]*rA(:,i);
        tAe(1:3,i) = rA(:,i)*cos(tma(1))+pA(:,i)*(dot(pA(:,i),rA(:,i)))*(1-cos(tma(1)))+...
            cross(rA(:,i),pA(:,i))*sin(tma(1));
        tA(1:3,i) = tAe(:,i)*cos(tma(2))+yA(:,i)*(dot(yA(:,i),tAe(:,i)))*(1-cos(tma(2)))+...
            cross(tAe(:,i),yA(:,i))*sin(tma(2));
%         tA(2,i) = -tA(2,i);
%         tA(3,i) = sqrt(1-tA(2,i)^2-tA(1,i)^2);
%         T(i,:) = (Tt(i) + (P(1)-P(i))*(pi*(dnoz/2/12)^2))*transpose(rA(:,i));
%         T(i,:) = (Tt(i) + (P(1)-P(i))*(pi*(dnoz/2/12)^2))*sqrt(1-tA(2,i)^2-tA(1,i)^2)*...
%             transpose(tA(:,i));
        T(i,:) = (1*Tt(i) + (P(1)-P(i))*(pi*(dnoz/2/12)^2))*transpose(tA(:,i));
%         yoyo = [cos(tma(2)),-sin(tma(2)),0;sin(tma(2)),cos(tma(2)),0;0,0,1]*...
%             %[cos(tma(1)),0,sin(tma(1));0,1,0;-sin(tma(1)),0,cos(tma(1))]*...
%             [1,0,0;0,cos(tma(1)),-sin(tma(1));0,sin(tma(1)),cos(tma(1))];
    end
%     Tx(i) = -Tt(i)*sin(p(i));
%     Tz(i) = Tt(i)*cos(p(i));
    Tr(i,1:3) = cross(rA(:,i),cross(rA(:,i),T(i,:)));
    Tr(i,1) = -Tr(i,1);
%     Tr(i,1:2) = -Tr(i,1:2);
%     T(i,:) = T(i,:)+Tr(i,:);

    %     wfr(i) = norm(T(i,:))/Isp;
    if t>=tb
        w(i) = dw;
        Tt(i) = 0;
    end
    
    %% Calculate Center of Mass and Stability
    
    wfr(i) = Tt(i)/Isp;
    w(i+1) = w(i) - wfr(i)*dt;
    wfrf(i) = wfr(i)/(OF+1);
    wf(i+1) = wf(i) - wfrf(i)*dt;
    wfro(i) = wfrf(i)*OF;
    wo(i+1) = wo(i) - wfro(i)*dt;
    lf(i+1) = (wf(i+1)/g)/(pi*rtank^2*df);
    lo(i+1) = (wo(i+1)/g)/(pi*rtank^2*dox);
    
    xfuel = xfuel0 + lfi - lf(i);
    xox = xox0 + loi - lo(i);
    cgf(i+1) = xfuel + lf(i)/2;
    cgo(i+1) = xox + lo(i)/2;
    if t<tb
        cg(i) = (cgd*dw + cgf(i)*wf(i) + cgo(i)*wo(i))/(w(i));
    else
        cg(i) = cgd;
        cg(i+1) = cgd; %To cover blind spots
    end
    if v(i,3) >=0 && dist(i)>RL
        stab(i) = (CP8(i)-cg(i))/(d/12);
    end
    
    %% Find velocities relative to wind and angle of attack of rocket
    
    vrwi(i,:) = vrw(i,:);
    if norm(omega(i,:)) == 0
        vrot(i,:) = [0,0,0];
    else
        vrot(i,:) = abs(CP8(i)-cg(i))*sin(acos(dot(rA(:,i),omega(i,:)/...
            norm(omega(i,:)))))*transpose(cross(rA(:,i),omega(i,:)))/...
            norm(cross(rA(:,i),omega(i,:)/norm(omega(i,:))));
    end

    vrw(i,:) = vrw(i,:)+vrot(i,:);
    aoa(i) = acos(dot((vrw(i,:)/norm(vrw(i,:))),rA(:,i)));
    if v(i,3)<0
        aoa(i) = 0;
    end
    vhat(i,:) = vrw(i,:)/norm(vrw(i,:));
    
%     if norm(omega(i,:)) == 0
%         vrotb(i,:) = [0,0,0];
%         vrotf(i,:) = [0,0,0];
%     else
%         vrotb(i,:) = vrw(i,:)+vrot(i,:)/abs(CP8(i)-cg(i))*(cgi(2)-cg(i));
%         vrotf(i,:) = vrw(i,:)+vrot(i,:)/abs(CP8(i)-cg(i))*(cgi(4)-cg(i));
        vrotb(i,:) = vrw(i,:)+vrot(i,:)/abs(CP8(i)-cg(i))*(CPi(i,2)-cg(i));
        vrotf(i,:) = vrw(i,:)+vrot(i,:)/abs(CP8(i)-cg(i))*(CPi(i,4)-cg(i));
%     end
    
%     vrot(i,:) = (cgi(:)-cg(i))*omega(i);
%     vrotx(i,:) = vrot(i,:)*cos(p(i));
%     vrotz(i,:) = vrot(i,:)*sin(p(i));
%     vrwrot(i,:) = ((vrotx(i,:)+vrw(i,1)).^2+(vrotz(i,:)+vrw(i,2)).^2).^0.5;
    
    %% 0 aoa Drag Coefficient Model
    
    mach(i) = norm(vrw(i,:))/vs(i);
    if i==1
        CDb6 = 0;
    end
    vrwrot(i,:) = [0,norm(vrotb(i,:)),0,norm(vrotf(i,:))];
    [CD(i),CDf(i),CDb(i),CDt(i),CDs(i),CDfb(i),CDff(i),CDfpro(i),CDll(i),CDb6] = ...
    getCD0(temp(i),mu0,temp0,rho(i),vrwrot(i,:),vrw(i,:),mach(i),Lc,tr,ft,Cr,Ct,S,d,Rl,n,...
    Afe2,machD,machF,dCDmax,aPro,Apro,Spro,aLL,All,Sll,Kb,dr,nb,t,tb,CDb6);
    Clpro(i) = 0.2; %% For now from CFD average
    
    %% Induced Drag Coefficient Model
    
    [CDab,CDaf] = getCDi(aoa(i),aoadd,aoadn,Rl,Ln,Afp,dfa,kfb,kbf,Afe,dist(i),RL,d);
    CDabs(i) = CDab;
    CDafs(i) = CDaf;
    
    %% Calculate Lift Coefficient
    
    cl8(i) = sum(CNa(i,:))*aoa(i)+Clpro(i);

    %% Get Total Drag Coefficient
    
%     CDa(i) = 5*sin(abs(aoa(i)));
    CDa(i) = 0;
%     CDa(i) = CDab+CDaf;
    cd8(i) = CD(i)+CDa(i);
    cd(i) = cd8(i);

    if cd8(i) < 0.3
        cd8(i) = 0.3;
    elseif cd8(i) > 1.6
        cd8(i) = 1.6;
    end
    
    %% Rolling Model
    
    Xf = CPi(i,4);

    rf = (r0/12+DR)/2+S/3;
    
    for fi=1:1:4
       PiB(fi,:) = [rf*sin(-pi/2+fi*pi/2),rf*cos(-pi/2+fi*pi/2),-Xf];
       Pi(fi,:) = transpose(RM*transpose(PiB(fi,:)))+x(i,:);
       Si(fi,:) = x(i,:)-Pi(fi,:);
%        Si(fi,:) = transpose(-RM*transpose(PiB(fi,:)));
       if norm(omega(i,:))~=0
           Vp(fi,:) = norm(omega(i,:))*norm(Si(fi,:))*sin(acos(dot((Si(fi,:)/norm(Si(fi,:))),...
               (omega(i,:)/norm(omega(i,:))))))*cross((Si(fi,:)/norm(Si(fi,:))),...
               (omega(i,:)/norm(omega(i,:))))/norm(cross((Si(fi,:)/norm(Si(fi,:))),...
               (omega(i,:)/norm(omega(i,:)))))+vrw(i,:);
       else
           Vp(fi,:) = vrw(i,:);
       end
%        li(fi,:) = [cos(-pi/2+fi*pi/2),sin(-pi/2+fi*pi/2),0];
       if fi==1
           li(fi,:) = transpose(pA(:,i));
       elseif fi==2
           li(fi,:) = transpose(yA(:,i));
       elseif fi==3
           li(fi,:) = transpose(-pA(:,i));
       elseif fi==4
           li(fi,:) = transpose(-yA(:,i));
       end
       Qc(1,i) = cos(cant/2);
       Qc(2:4,i) = sin(cant/2)*cross(li(fi,:),rA(:,i))/norm(cross(li(fi,:),rA(:,i)));
       Rc = [1-2*Qc(3,i)^2-2*Qc(4,i)^2, 2*Qc(2,i)*Qc(3,i)-2*Qc(1,i)*Qc(4,i),...
           2*Qc(2,i)*Qc(4,i)+2*Qc(1,i)*Qc(3,i);...
           2*Qc(2,i)*Qc(3,i)+2*Qc(1,i)*Qc(4,i),...
           1-2*Qc(2,i)^2-2*Qc(4,i)^2, 2*Qc(3,i)*Qc(4,i)-2*Qc(1,i)*Qc(2,i);...
           2*Qc(2,i)*Qc(4,i)-2*Qc(1,i)*Qc(3,i),...
           2*Qc(3,i)*Qc(4,i)+2*Qc(1,i)*Qc(2,i), 1-2*Qc(2,i)^2-2*Qc(3,i)^2];
       aoafi(i,fi) = pi/2-acos(dot((Vp(fi,:)/norm(Vp(fi,:))),Rc*transpose(li(fi,:))));
       Vpn(i,fi) = norm(Vp(fi,:));
    end
    aoaf(i) = sum(aoafi(i,:))/n;
    Vf(i) = sum(Vpn(i,:))/n;
%     CR(i) = aoaf(i)/750;
    LFr(i) = (CNaf/n)*aoaf(i);

    rf = (r0/12+DR)/2;
    CNa0 = 2*pi/beta;
    Clf(i) = n*(MAC+rf)*(CNaf/n)*cant/(d/12);
    term = (Cr+Ct)/2*rf^2*S+(Cr+2*Ct)/3*rf*S^2+(Cr+3*Ct)/12*S^3;
    if mach(i)<1
        rollrate(i) = Aref*beta*norm(vrw(i,:))*MAC*(CNaf/n)*cant/(2*pi*term);
    elseif mach(i)>=1
        rollrate(i) = Aref*beta*norm(vrw(i,:))*MAC*(CNaf/n)*cant/(2*pi*term);
    end
    Cdf(i) = n*CNa0*rollrate(i)*term/Aref/(d/12)/norm(vrw(i,:));
    CR(i) = Clf(i)-Cdf(i);
    rf = (r0/12+DR)/2+S/3;
    
    %% Calculate Drag Force
    
    %L(i,:) = cl8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*[cos(p(i)),sin(p(i))];
    
%     Lx(i) = L(i)*cos(p(i));
%     Lz(i) = L(i)*sin(p(i));
    
%     D(i,:) = cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*[sin(phi(i)),-cos(phi(i))];
% %     Dx(i) = D(i)*sin(phi(i));
% %     Dz(i) = -D(i)*cos(phi(1));
    if v(i,3) > 0
%         D(i,:) = cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*[sin(phi(i)),-cos(phi(i))];
        D(i,:) = -cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*rA(:,i);
%         D(i,:) = -cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*vhat(i,:);
%         if v(i,2) == 0 && t>tb
%             tff = t;
%             itf = i;
%         end
    elseif norm(v(i,:))==0
        D(i,:) = [0,0,0];
    else
        phi(i) = -pi-phi(i);
        %D(i) = cd(i)/2*Aref*rho(i)*(vrw(i))^2;
%         D(i,:) = cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*[sin(phi(i)),-cos(phi(i))];
        D(i,:) = -cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*vhat(i,:); 
        %         Dx(i) = D(i)*sin(phi(i));
        %         Dz(i) = -D(i)*cos(phi(i));
        if v(i-1,2) >=0 && ejection == 0
            Iapogee2 = i;
        end
        if recovery == 0 && dist(i)>RL
            tff = t;
            itf = i;
%             disp(i+1)
            break
        end
        if zs(i) > 1500
            cd8(i) = 1.6;
            d = 24;
            Aref = pi*(d/12/2)^2;
%             D(i,:) = cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*[sin(phi(i)),-cos(phi(i))];
            D(i,:) = -cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*vhat(i,:);
            %             D(i) = cd(i)/2*A*rho(i)*(vrw(i))^2;
            %             Dx(i) = D(i)*sin(phi(i));
            %             Dz(i) = -D(i)*cos(phi(i));
        else
%             if cd8(i-1) == 1.6
%                 Reac(i,2) = 800;
%                 ejection = 0;
%             end
            cd8(i) = 2.2;
            d = 60;
            Aref = pi*(d/12/2)^2;
%             D(i,:) = cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*[sin(phi(i)),-cos(phi(i))];
            D(i,:) = -cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*vhat(i,:);
            %             D(i) = cd(i)/2*A*rho(i)*(vrw(i))^2;
            %             Dx(i) = D(i)*sin(phi(i));
            %             Dz(i) = -D(i)*cos(phi(i));
        end
    end
    
    %% Calculate Lift Force
    
    L(i,:) = cl8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*...
        cross(rA(:,i),cross(rA(:,i),vrw(i,:)/norm(vrw(i,:))))/...
        norm(cross(rA(:,i),cross(rA(:,i),vrw(i,:)/norm(vrw(i,:)))));
    if v(i,3) < 0
        L(i,:) = [0,0,0];
    end    
    
    %% Calculate Weight
    
    if t<=tb
        W(i,:) = [0,0,-w(i)];
    else
        W(i,:) = [0,0,-dw];
    end
    
    %% Calculate Fictitious Forces
    
    Eomega1 = [0,Eomega,0];
    C(i,:) = 2*m*cross(Eomega1,v(i,:))-m*cross(Eomega1,cross(Eomega1,x(i,:))); %% Correct this
    %C(i,:) = [0,0];
    
    %% Calculate Rail Friction
    
    if s(i) > 0 && dist(i) < RL && t<tb
%         Frfk(i) = muk*w(i)*sin(abs(p(i)))+2*muk/Lll*((norm(T(i,:))-w(i)*cos(p(i))-norm(D(i)))*...
%             (r0/12+rll)+w(i)*sin(abs(p(i)))*(aLL/12-cg(i)));
        Frfk(i) = muk*(w(i)*sin(abs(p(i)))+norm(Tr(i,:)));
        Frk(i,:) = abs(Frfk(i))*[sin(p(i)),0,-cos(p(i))];
    elseif dist(i) <= 0 && t<tb
%         Frfs(i) = mus*w(i)*sin(abs(p(i)))+2*mus/Lll*((norm(T(i,:))-w(i)*cos(p(i)))*...
%             (r0/12+rll)+w(i)*sin(abs(p(i)))*(aLL/12-cg(i)));
        Frfs(i) = mus*(w(i)*sin(abs(p(i)))+norm(Tr(i,:)));
        Frs(i,:) = abs(Frfs(i))*[sin(p(i)),0,-cos(p(i))];
%         Frfk(i) = muk*w(i)*sin(abs(p(i)))+2*muk/Lll*((norm(T(i,:))-w(i)*cos(p(i))-norm(D(i)))*...
%             (r0/12+rll)+w(i)*sin(abs(p(i)))*(aLL/12-cg(i)));
        Frfk(i) = muk*(w(i)*sin(abs(p(i)))+norm(Tr(i,:)));
        Frk(i,:) = abs(Frfk(i))*[sin(p(i)),0,-cos(p(i))];
    end
    
    %% Calculate Reaction Forces and Moments due to Recovery Shock
    
    if Iapogee2 ~= 0 && ejection==0 && cd8(i) ~= 2.2
        Shock = [1700/3,1700];
%         Shock = [0,0];
        Reac(i,:) = Shock(1)*[-sign(v(i,1))*cos(p(i)),0,-sin(p(i))]; %Make perpendicular to axis
        Reac(i,:) = Reac(i,:)+Shock(2)*[-sin(p(i)),0,cos(p(i))]; %Axial
        Mreac = [0,-sign(p(i))*Shock(1)*(cg(i)-47/12),0];
        ejection = 1;
        % [-sin(p(i)),cos(p(i))]
    else
        Mreac = [0,0,0];
    end
    
    %% Calculate Total Force on Rocket
    
    if dist(i) > RL || v(i,3) < 0 % dist(i) > RL-Rl+cg(i) or dist(i) > RL for OTRS?
        F(i,:) = T(i,:)+W(i,:)+D(i,:)+L(i,:)-C(i,:)+Reac(i,:);
        %         Fx(i) = Tx(i)+Dx(i)+Lx(i);
        %         Fz(i) = Tz(i)+Dz(i)+Lz(i)-w(i);
        if v(i,3)>-1
%             AF(i) = norm(T(i,:))-norm(D(i,:))-norm(W(i,:))*cos(p(i));
%             AF(i) = dot(F(i,:),[-sin(p(i)),cos(p(i))]);
%             NF(i) = -sign(aoa(i))*norm(L(i,:))+norm(W(i,:))*sin(p(i));
%             NF2(i,:) = cross([F(i,1),0,F(i,2)],[-sin(p(i)),0,cos(p(i))]);
%             NF(i) = dot(F(i,:),[-cos(p(i)),-sin(p(i))]);
        end
%         if Reac(i,2)==800
%             I_main = i;
%             AF(i) = norm(W(i,:))-norm(D(i,:))-norm(Reac(i,:));
%             NF(i) = 0;
%         end
%     elseif ((norm(T(i,:))-w(i)*cos(p(i)))>abs(Frfs(end)) || s(i) > 0) && t<tb  %Moving on rail
    elseif ((norm(T(i,:))-w(i)*cos(p(i))-norm(D(i,:)))>abs(Frfs(end)) || s(i) > 0) && t<tb
%         Reac(i,:) = [norm(W(i,:))*sin(p(i))*cos(p(i)),0,norm(W(i,:))*sin(p(i))^2];
        F(i,:) = dot((T(i,:)+W(i,:)+D(i,:)+Frk(i,:)),rA(:,i))*transpose(rA(:,i));
%         Fx(i) = Tx(i)+Dx(i)+Lx(i);
%         Fz(i) = Tz(i)+Dz(i)+Lz(i)-w(i);
%         AF(i) = norm(T(i,:))-norm(D(i,:))-norm(W(i,:))*cos(p(i));
%         AF2(i) = dot(F(i,:),[-sin(p(i)),cos(p(i))]);
%         NF(i) = norm(W(i,:))*sin(p(i))-sign(p(i))*norm(Reac(i,:));
        OTRs = s(i);
        OTRi = i;
%     elseif ((norm(T(i,:))-w(i)*cos(p(i)))<=abs(Frfs(i)) && dist(i) <= 0) && t<tb
    elseif ((norm(T(i,:))-w(i)*cos(p(i)))<=abs(Frfs(i)) && dist(i) <= 0) && t<tb
%         Reac(i,:) = -(T(i,:)+W(i,:)+Frs(i,:));
%         F(i,:) = T(i,:)+W(i,:)+Reac(i,:)+Frs(i,:);
        F(i,:) = [0,0,0];
%         Fx(i) = 0;
%         Fz(i) = 0;
        L(i,:) = [0,0,0];
    end
    
    %% Calculate Linear Motion Variables
    
%     ax(i) = Fx(i)/m;
%     az(i) = Fz(i)/m;
    a(i,:) = F(i,:)/m;
%     vx(i+1) = vx(i) + ax(i)*dt;
%     vz(i+1) = vz(i) + az(i)*dt;
    v(i+1,:) = v(i,:) + a(i,:)*dt;
    x(i+1,:) = x(i,:) + v(i,:)*dt;
%     z(i+1) = z(i) + vz(i)*dt;
    dist(i+1) = norm(x(i+1,:)-x(i,:))+dist(i);
    s(i+1) = norm(v(i+1,:));
    aA(i) = AF(i)/m;
    aN(i) = NF(i)/m;
    
    %% Check for Instability

    if CP8(i)<cg(i) && dist(i)>RL-Rl+cg(i) && v(i,3)>0 && abs(aoa(i))<15*pi/180
%         disp(i);
        %error('Unstable');
    end
    
    %% Calculate Moments of Inertia
    
    If(i) = (1/12)*(wf(i))*(3*rtank^2 + lf(i)^2) + (wf(i))*(cg(i)-cgf(i))^2;
    Io(i) = (1/12)*(wo(i))*(3*rtank^2 + lo(i)^2) + (wo(i))*(cg(i)-cgo(i))^2;
    if t>tb
        Ixx(i) = Id;
        Iyy(i) = Id;
        Izz(i) = Idz;
    else
        Ixx(i) = Id + (dw)*(cg(i)-cgd)^2 + If(i) + Io(i);
        Iyy(i) = Ixx(i);
        Izz(i) = Idz+(1/2)*((wf(i))*rtank^2+(wo(i))*rtank^2);
    end
    I(i,:) = [Ixx(i),Iyy(i),Izz(i)];
    I0 = diag(I(i,:),0);
    I0 = RM*I0*transpose(RM);
%     I0 = [yA(:,i),pA(:,i),rA(:,i)]*I0*transpose([yA(:,i),pA(:,i),rA(:,i)]);
%     I0 = transpose(RM)*I0*RM;
%     I0 = transpose(RM)*I0*RM;
    I0Rec{i} = I0;
    if i ~= 1
        Idot = (I0Rec{i}-I0Rec{i-1})/dt;
    else
        Idot = zeros(3,3);
    end
    
    %% Calculate Moments 
    
    cgp(i) = (cgf(i)*wf(i)+cgo(i)*wo(i))/(wf(i)+wo(i));
    Mt(i,:) = -(wfr(i)/g)*((Rl-cg(i))^2-(cgp(i)-cg(i))^2)*... % Causes problems w/ wind
        transpose(RM*diag([1 1 0],0)/RM*transpose(omega(i,:))); % directions and stuff
    Mt(i)=0;
%     LM(i,:) = [cg(i)-CPi(i,1),cg(i)-CPi(i,2),CPi(i,3)-cg(i),CPi(i,4)-cg(i)];
%     MN(i) = Q(i)*Aref*sum(CNa(i,:).*aoai(i,:).*LM(i,:));
    MN(i,:) = -abs(CP8(i) - cg(i))*norm(L(i,:))*cross(rA(:,i),vrw(i,:)/norm(vrw(i,:)))/...
        norm(cross(rA(:,i),vrw(i,:)/norm(vrw(i,:))));
%     Cm(i) = -stab(i)*cl8(i);
%     MNcm(i) = Cm(i)*0.5*rho(i)*norm(vrw(i,:))^2*Aref*(d/12);
    %summ = CNa(i,1).*(CPi(1)-cg(i)).^2+CNa(i,2).*(CPi(2)-cg(i)).^2+...
    %      CNa(i,3).*(CPi(3)-cg(i)).^2+CNa(i,4).*(CPi(4)-cg(i)).^2;
%     thing = CNa(i,:).*(CPi(i,:)-cg(i)).^2;
%     summ = trace(CNa(i,:).*(CPi(i,:)-cg(i)).^2);
%     summ = sum(CNa(i,:).*(CPi(i,:)-cg(i)).^2);
    %summ = CNa(i,:).*(CPi(i,:)-cg(i)).^2;
%     Md(i) = 0.5*rho(i)*norm(vrw(i,:))*Aref*omega(i)*summ;
%     Md(i) = 0;
    Mll(i,1:3) = -0.5*rho(i)*(norm(vrw(i,:)))^2*...
        (CDll(i))*Aref*(r0/12)*transpose(pA(:,i)); %+-Depends on which side lug is on
    Mpro(i,1:3) = 0.5*rho(i)*(norm(vrw(i,:)))^2*...
        (CDfpro(i))*Aref*(r0/12)*transpose(pA(:,i)); %+-Depends on protuberance side
    Dpro(i) = 0.5*rho(i)*(norm(vrw(i,:)))^2*(CDfpro(i))*Aref*(r0/12);
    MT(i,:) = -(Rl-cg(i))*cross(rA(:,i),Tr(i,:));
%     MT(i,:) = [0,0,0];
%     M(i) = MN(i)+Mt(i)+Md(i)+Mll(i)+Mpro(i)+Mreac;
%     Mr(i,:) = CR(i)/2*rho(i)*Aref*(rollrate(i)*rf)^2*(d/12)*transpose(rA(:,i));
    Mr(i,:) = -CR(i)*LFr(i)*rf*transpose(rA(:,i));
%     Mr(i,:) = [0,0,0];
    M(i,:) = MN(i,:)+Mt(i,:)+Mll(i,:)+Mpro(i,:)+Mreac+MT(i,:)+Mr(i,:);
%     alpha(i,:)=(I0)\transpose(M(i,:)-transpose(Idot*transpose(omega(i,:))));

    %% Calculate Rotational Motion Variables

    alpha(i,:)=(I0)\(transpose(M(i,:))-...
        transpose(cross(omega(i,:),I0*transpose(omega(i,:)))));
    if dist(i)<RL || v(i,3)<0
        alpha(i,:) = [0,0,0];
    end
    omega(i+1,:) = omega(i,:) + alpha(i,:)*dt;
    y(i+1) = y(i) + omega(i,1)*dt;
    p(i+1) = p(i) + omega(i,2)*dt;
    r(i+1) = r(i) + omega(i,3)*dt;
    
    %% Update Quaternion-related Values
    
    q(:,i+1) = [cos(r(i+1)/2)*cos(p(i+1)/2)*cos(y(i+1))+sin(r(i+1)/2)*sin(p(i+1)/2)*sin(y(i+1));...
     sin(r(i+1)/2)*cos(p(i+1)/2)*cos(y(i+1))-cos(r(i+1)/2)*sin(p(i+1)/2)*sin(y(i+1));...
     cos(r(i+1)/2)*sin(p(i+1)/2)*cos(y(i+1))+sin(r(i+1)/2)*cos(p(i+1)/2)*sin(y(i+1));...
     cos(r(i+1)/2)*cos(p(i+1)/2)*sin(y(i+1))-sin(r(i+1)/2)*sin(p(i+1)/2)*cos(y(i+1))];
%     q(:,i+1) = q(:,1)+qdot(:,i)*dt;
    qdot(:,i+1) = 0.5*([0,-omega(i+1,1),-omega(i+1,2),-omega(i+1,3);...
        omega(i+1,1),0,omega(i+1,3),-omega(i+1,2);...
        omega(i+1,2),-omega(i+1,3),0,omega(i+1,1);...
        omega(i+1,3),omega(i+1,2),-omega(i+1,1),0]*q(:,i+1));
%     r(i+1) = atan2(2*(q(1,i+1)*q(2,i+1)+q(3,i+1)*q(4,i+1)),1-2*(q(2,i+1)^2+q(3,i+1)^2));
%     p(i+1) = asin(2*(q(1,i+1)*q(3,i+1)-q(4,i+1)*q(2,i+1)));
%     y(i+1) = atan2(2*(q(1,i+1)*q(4,i+1)+q(2,i+1)*q(3,i+1)),1-2*(q(3,i+1)^2+q(4,i+1)^2));
    RM = [1-2*q(3,i+1)^2-2*q(4,i+1)^2, 2*q(2,i+1)*q(3,i+1)-2*q(1,i+1)*q(4,i+1),...
        2*q(2,i+1)*q(4,i+1)+2*q(1,i+1)*q(3,i+1);...
        2*q(2,i+1)*q(3,i+1)+2*q(1,i+1)*q(4,i+1),...
        1-2*q(2,i+1)^2-2*q(4,i+1)^2, 2*q(3,i+1)*q(4,i+1)-2*q(1,i+1)*q(2,i+1);...
        2*q(2,i+1)*q(4,i+1)-2*q(1,i+1)*q(3,i+1),...
        2*q(3,i+1)*q(4,i+1)+2*q(1,i+1)*q(2,i+1), 1-2*q(2,i+1)^2-2*q(3,i+1)^2];
    RM(1:2,:) = -RM(1:2,:);
    RMrec{i+1} = RM;
    yA(:,i+1) = RM*[1;0;0];
    pA(:,i+1) = RM*[0;1;0];
    rA(:,i+1) = RM*[0;0;1];
%     rA(1,i+1) = -rA(1,i+1);
%     rA(2,i+1) = -rA(2,i+1);
    
    %% Miscellaneous Values
    
    gs(i) = g;
    ms(i) = m;
%     if dist(i) > RL || v(i,2) < 0 && v(i,2)>0
%         NF2(i) = sign(alpha(i))*norm(cross(F(i,:),[-sin(p(i)),cos(p(i))]));
%     end

    Q(i) = 0.5*rho(i)*(norm(vrw(i,:)))^2;
    Qa(i) = Q(i)*aoa(i);
    QCNa(i) = Q(i)*sum(abs(CNa(i,:)));
    
    %% Iterate
    
    iters(i) = i;
    i = i + 1;
    
end    

%% Results
% Find apogee, max velocity, max acceleration
% as well as the times at which they happen

if recovery ~= 0
    range_bound_right = max(x(:,1));
    range_bound_left = min(x(:,1));
    landing_position = x(itf,1);
end
[apogee,Iapogee] = max(zs);
% fprintf("Right bound: %f\nLeft bound: %f\nLanding position: %f\n",...
%     range_bound_right,range_bound_left,landing_position);

[AF_max,IAF_max] = max(abs(AF));
[NF_max,INF_max] = max(abs(NF));
%fprintf("Max A = %f at t = %.2f\nMax N = %f at t = %.2f\n",AF_max,IAF_max*dt,NF_max,INF_max*dt);
T_mag = zeros(t_steps,1);
for it = 1:1:t_steps
    T_mag(it) = norm(T(it,:));
    L_mag(it) = norm(L(it,:));
    D_mag(it) = norm(D(it,:));
end
T_max = max(T_mag);
[D_max,I_Dmax] = max(D_mag);
L_max = max(L_mag);
[M_max,I_Mmax] = max(abs(M));
[maxQ,I_maxQ] = max(Q);
[maxQa,I_maxQa] = max(abs(Qa));
[maxQCNa,I_maxQCNa] = max(abs(Qa));
FinLift = 0.5*rho(I_maxQCNa)*CNa(I_maxQCNa,4)*aoa(I_maxQCNa)*Aref*norm(vrw(I_maxQCNa,:))^2;
% I_maxQCNa = Iapogee2;
% I_maxQCNa = 2;
% I_maxQCNa = I_main;
OTR_Tavg = sum(T(1:OTRi,3))/OTRi;
% OTR_Tavg = sum(Tt(1:floor(tb*100)))/floor(tb*100);
deployment_v = s(Iapogee);
% fprintf("Apogee: %f\nDeployment Velocity: %f\nOTR Speed: %f\n",apogee,deployment_v,OTRs);
% [ymax,t1I] = max(ys);
% [vmax,t2I] = max(vs);
% [amax,t3I] = max(as);
% t1 = t1I*dt;
% t2 = t2I*dt;
% t3 = t3I*dt;


m1 = 9.62; %Nose
m2 = 28.48+2.5;
m3 = 8.87+wo(1);
m4 = 2.26;
m5 = 8.12+wf(1)+2.5;
m6 = 9.35; %Right Above Thrust Bulkhead
m7 = 19.39;
masses = [m1,m2,m3,m4,m5,m6,m7]/gs(I_maxQCNa);
distances = [0,47,36.43,loi*12,5.5,lfi*12,15.5,19.3,0];
RB = distances(2);
% DfbT = norm(D(I_maxQCNa,:))*(CDfb(I_maxQCNa)+CDfpro(I_maxQCNa)+CDll(I_maxQCNa)+...
%     CDabs(I_maxQCNa)/betas(I_maxQCNa))/cd8(I_maxQCNa); %Model protuberances separately later
DfbT = norm(D(I_maxQCNa,:))*0.6; %Model protuberances separately later
% Dst = norm(D(I_maxQCNa,:))*(CDt(I_maxQCNa)+CDs(I_maxQCNa))/cd8(I_maxQCNa);
Dst = norm(D(I_maxQCNa,:))*0.1;
Reac1(I_maxQCNa,2) = dot(Reac(I_maxQCNa,:),[-sin(p(I_maxQCNa)),0,cos(p(I_maxQCNa))]);
Reac1(I_maxQCNa,1) = dot(Reac(I_maxQCNa,:),[-cos(p(I_maxQCNa)),0,-sin(p(I_maxQCNa))]);

% for mi = 1:1:length(masses)
% %     Dfb(mi) = DfbT*sum(distances(1:mi))/sum(distances)+Dst;
%     Dfb(mi) = DfbT+Dst-Reac1(I_maxQCNa,2); %Recovery
%     if mi == length(masses)
% %         Dbff = norm(D(I_maxQCNa,:))*(CDff(I_maxQCNa)+CDb(I_maxQCNa)+...
% %             CDafs(I_maxQCNa)/betas(I_maxQCNa))/cd8(I_maxQCNa);
%         Dbff = norm(D(I_maxQCNa,:))*0.3;
% %         Dfb(mi) = Dfb(mi)+Dbff-norm(T(I_maxQCNa,:));
%         Dfb(mi) = Dfb(mi)+Dbff-norm(T(I_maxQCNa,:)); %Recovery
%     end
%     if mi >= 4
%         weight = w(I_maxQCNa)*cos(p(I_maxQCNa));
%     else
%         weight = 0;
%     end
%     FAL(mi+1) = sum(masses(1:mi))*abs(aA(I_maxQCNa))+Dfb(mi)+sum(masses(1:mi))*...
%         gs(I_maxQCNa)*cos(p(I_maxQCNa));
% %     FAL(mi+1) = sum(masses(1:mi))*abs(aA(I_maxQCNa))+Dfb(mi);%+weight;
% end
FAL2(length(masses)+1) = 0;
for mi = length(masses):-1:1
%     Dfb(mi) = DfbT*sum(distances(1:mi))/sum(distances)+Dst;
%     Dfb(mi) = DfbT+Dst;
    if mi == 1
        Di(mi) = DfbT+Dst;
%         Dfb2(mi) = 0;
    elseif mi == length(masses)
        Di(mi) = norm(D(I_maxQCNa,:))*0.3;
    else
        Di(mi) = 0;
    end
    if mi == 1
        Reaci(mi) = Reac1(I_maxQCNa,2);
    else
        Reaci(mi) = 0;
    end
    if mi == length(masses)
%         Ti(mi) = norm(T(I_maxQCNa-1,:)); %Max q-alpha
        Ti(mi) = norm(T(I_maxQCNa,:));
    else
        Ti(mi) = 0;
    end
%     Dbff2 = norm(D(I_maxQCNa,:))*(CDff(I_maxQCNa)+CDb(I_maxQCNa)+...
%         CDafs(I_maxQCNa)/betas(I_maxQCNa))/cd8(I_maxQCNa);
    FIg(mi) = masses(mi)*(aA(I_maxQCNa)+gs(I_maxQCNa));
    FAL2(mi) = FAL2(mi+1)+Ti(mi)-Di(mi)-FIg(mi)+Reaci(mi);
%     FAL2(mi) = FAL2(mi+1)-Ti(mi)+Di(mi)-FIg(mi)-Reaci(mi);
%     FAL2(mi) = -sum(masses(mi:end))*abs(aA(I_maxQCNa))-Dfb22(mi);%+weight;
end
% FAL2(length(masses)+2) = 0;

% d = 7;
% Aref = pi*(d/12/2)^2;
% Lnose(1) = 0;
% Lbody(1) = 0;
% Ltail(1) = 0;
% Lfins(1) = 0;
% cgs(length(masses)+1,1) = 0;
% % BM(length(masses)+1) = 0;
% aN2 = ((-0.5*aoa(I_maxQCNa)*sum(CNa(I_maxQCNa,:))*Aref*rho(I_maxQCNa)*norm(vrw(I_maxQCNa,:))^2+...
%     Reac1(I_maxQCNa,1))/(w(I_maxQCNa)/gs(I_maxQCNa)));
% FVL(1) = 0;
% BM(1) = 0;
% Reac2 = 0;
% sum1 = zeros(1,length(masses));
% sum2 = zeros(1,length(masses));
% sum3 = zeros(1,length(masses));
% for mi = 1:1:length(masses)
%     cgs(mi,1) = sum(distances(1:mi))+distances(mi+1)/2;
% %     cgs(mi) = Rl*12-cgs(mi);
%     if mi == length(masses)
%         Ltail(mi) = CNa(I_maxQCNa,3)*aoai(I_maxQCNa,3)/2*Aref*rho(I_maxQCNa)*vrwrot(I_maxQCNa,3)^2;
%         Lfins(mi) = CNa(I_maxQCNa,4)*aoai(I_maxQCNa,4)/2*Aref*rho(I_maxQCNa)*vrwrot(I_maxQCNa,4)^2;
%     else
%         Ltail(mi) = 0;
%         Lfins(mi) = 0;
%     end
%     if mi>=3
%         Lbody(mi) = CNa(I_maxQCNa,2)*aoai(I_maxQCNa,2)/2*Aref*rho(I_maxQCNa)*...
%             vrwrot(I_maxQCNa,2)^2;
%     else
%         Lbody(mi) = 0;
%     end
%     if mi>=1
%         Lnose(mi) = CNa(I_maxQCNa,1)*aoai(I_maxQCNa,1)/2*Aref*rho(I_maxQCNa)*...
%             vrwrot(I_maxQCNa,1)^2;
%         Reac2(mi) = Reac1(I_maxQCNa,1);
%     end
% %     sum1(1,mi) = 0;
% %     sum2(1,mi) = 0;
% %     sum3(1,mi) = 0;
%     for mis = 1:1:mi
%         sum1(mis+1,mi) = sum1(mis,mi)+masses(mis)*(sum(distances(1:mi+1))-cgs(mis,1));
%         sum2(mis+1,mi) = sum2(mis,mi)+masses(mis)*(sum(distances(1:mi+1))-cgs(mis,1))*...
%             ((cg(I_maxQCNa))*12-cgs(mis,1));
%         sum3(mis+1,mi) = sum3(mis,mi)+masses(mis)*(cg(I_maxQCNa)-cgs(mis,1));
%     end
%     FVL(mi+1) = -Lnose(mi)-Lbody(mi)+Ltail(mi)-Lfins(mi)+aN2*sum(masses(mi))-...
%         alpha(I_maxQCNa)*sum3(mi+1,mi);
%     BM(mi+1) = Lnose(mi)/12*(sum(distances(1:mi+1))-12*(CPi(I_maxQCNa,1)))+...
%         Lbody(mi)/12*(sum(distances(1:mi+1))-12*CPi(I_maxQCNa,2))+...
%         Ltail(mi)/12*(sum(distances(1:mi+1))-12*CPi(I_maxQCNa,3))+...
%         Lfins(mi)/12*(sum(distances(1:mi+1))-12*CPi(I_maxQCNa,4))+...
%         Reac2(mi)/12*(sum(distances(1:mi+1))-12*Ln)+...
%         aN2*sum1(mi+1,mi)/12+alpha(I_maxQCNa)*sum2(mi+1,mi)/144;
% end
% BM(mi+2) = 0;
% 
% d = 7;
% Aref = pi*(d/12/2)^2;
% Lnose = 0;
% Lbody = 0;
% Ltail = CNa(I_maxQCNa,3)*aoai(I_maxQCNa,3)/2*Aref*rho(I_maxQCNa)*vrwrot(I_maxQCNa,3)^2;
% Lfins = CNa(I_maxQCNa,4)*aoai(I_maxQCNa,4)/2*Aref*rho(I_maxQCNa)*vrwrot(I_maxQCNa,4)^2;
% cgs(length(masses)+1,2) = 0;
% BM2(length(masses)+1) = 0;
% Reac2 = 0;
% for mi = length(masses):-1:1
%     cgs(mi,2) = sum(distances(mi+2:end))+distances(mi+1)/2;
% %     cgs(mi) = Rl*12-cgs(mi);
%     if mi<=3
%         Lbody = CNa(I_maxQCNa,2)*aoai(I_maxQCNa,2)/2*Aref*rho(I_maxQCNa)*...
%             vrwrot(I_maxQCNa,2)^2;
%     end
%     if mi==1
%         Lnose = CNa(I_maxQCNa,1)*aoai(I_maxQCNa,1)/2*Aref*rho(I_maxQCNa)*...
%             vrwrot(I_maxQCNa,1)^2;
%         Reac2 = Reac1(I_maxQCNa,1);
%     end
%     sum1 = 0;
%     sum2 = 0;
%     sum3 = 0;
%     for mis = length(masses):-1:mi
%         sum1 = sum1+masses(mis)*(sum(distances(mi+1:end))-cgs(mis,2))/12;
%         sum2 = sum2+masses(mis)*((sum(distances(mi+1:end))-cgs(mis,2))*...
%             ((Rl-cg(I_maxQCNa))*12-cgs(mis,2)))/144;
%         sum3 = sum3+masses(mis)*((Rl-cg(I_maxQCNa))*12-cgs(mis,2));
%     end
%     FVL(mi) = Lnose+Lbody+Ltail+Lfins-aN(I_maxQCNa)*sum(masses(length(masses):-1:mi))-...
%         alpha(I_maxQCNa)*sum(masses(length(masses):-1:mi)*((Rl-cg(I_maxQCNa))*12-...
%         cgs(length(masses):-1:mi,2))/12);
%     BM2(mi) = Lnose/12*(sum(distances(mi+1:end))-12*(Rl-CPi(I_maxQCNa,1)))+...
%         Lbody/12*(sum(distances(mi+1:end))-12*(Rl-CPi(I_maxQCNa,2)))+...
%         Ltail/12*(sum(distances(mi+1:end))-12*(Rl-CPi(I_maxQCNa,3)))+...
%         Lfins/12*(sum(distances(mi+1:end))-12*(Rl-CPi(I_maxQCNa,4)))+...
%         Reac2/12*(sum(distances(mi+1:end))-12*(Rl-Ln))-...
%         aN(I_maxQCNa)*sum1-alpha(I_maxQCNa)*sum2;
% end

d = 7;
Aref = pi*(d/12/2)^2;
Lnose(length(masses)) = 0;
Lbody(length(masses)) = 0;
Ltail(length(masses)) = -CNa(I_maxQCNa,3)*aoai(I_maxQCNa,3)/2*Aref*...
    rho(I_maxQCNa)*vrwrot(I_maxQCNa,3)^2;
Lfins(length(masses)) = -CNa(I_maxQCNa,4)*aoai(I_maxQCNa,4)/2*Aref*...
    rho(I_maxQCNa)*vrwrot(I_maxQCNa,4)^2;
cgs(length(masses)+1,2) = 0;
% BM2(length(masses)+1) = 0;
Reac2(length(masses)) = 0;
% aN2 = ((-0.5*aoa(I_maxQCNa)*sum(CNa(I_maxQCNa,:))*Aref*rho(I_maxQCNa)*norm(vrw(I_maxQCNa,:))^2+...
%     Reac1(I_maxQCNa,1))/(w(I_maxQCNa)/gs(I_maxQCNa)));
FVL(length(masses)+1) = 0;
BM3(length(masses)+1) = 0;

aN2 = aN(I_maxQCNa);
% distances(2) = 0;
% masses(1) = 0;
% aN2 = aN(I_maxQCNa)*(w(I_maxQCNa)/gs(I_maxQCNa))/sum(masses);

xCPi = zeros(length(masses),4);
for mi = length(masses):-1:1
    cgs(mi,2) = sum(distances(mi+2:end))+distances(mi+1)/2;
    xcgi(mi) = cgs(mi,2)/12-sum(distances(mi+2:end))/12;
    ls(mi) = distances(mi+1)/12;
    if mi<length(masses)
        Ltail(mi) = 0;
        Lfins(mi) = 0;
    elseif mi==length(masses)
        xCPi(mi,3) = Rl-CPi(I_maxQCNa,3)-sum(distances(mi+2:end))/12;
        xCPi(mi,4) = Rl-CPi(I_maxQCNa,4)-sum(distances(mi+2:end))/12;
        Lcp(mi) = Ltail(mi)*xCPi(mi,3)+Lfins(mi)*xCPi(mi,4);
    end
    if mi==4
        Lbody(mi) = -CNa(I_maxQCNa,2)*aoai(I_maxQCNa,2)/2*Aref*rho(I_maxQCNa)*...
            vrwrot(I_maxQCNa,2)^2;
        xCPi(mi,2) = Rl-CPi(I_maxQCNa,2)-sum(distances(mi+2:end))/12;
        Lcp(mi) = Lbody(mi)*xCPi(mi,2);
    else
        Lbody(mi) = 0;
    end
    FIi(mi) = masses(mi)*(-aN2+alpha(I_maxQCNa)*((Rl-cg(I_maxQCNa))-cgs(mi,2)/12));
    if mi==2 
        Reac2(mi) = Reac1(I_maxQCNa,1);
%         FIi(mi) = FIi(mi) -...
%             masses(mi-1)*(-aN2+alpha(I_maxQCNa)*((Rl-cg(I_maxQCNa))-cgs(mi-1,2)/12));
    else
        Reac2(mi) = 0;
    end
    if mi==1
        Lnose(mi) = -CNa(I_maxQCNa,1)*aoai(I_maxQCNa,1)/2*Aref*rho(I_maxQCNa)*...
            vrwrot(I_maxQCNa,1)^2;
        xCPi(mi,1) = Rl-CPi(I_maxQCNa,1)-sum(distances(mi+2:end))/12;
        Lcp(mi) = Lnose(mi)*xCPi(mi,1);
%         Lnose(mi) = 0;
%         Lcp(mi) = 0;
%         FIi(mi) = 0;
    else
        Lnose(mi) = 0;
    end
    Wi(mi) = masses(mi)*gs(I_maxQCNa)*sin(p(I_maxQCNa));
    Li(mi) = Lnose(mi)+Lbody(mi)+Ltail(mi)+Lfins(mi);
%     Li(mi) = 0;
%     Lcp(mi) = 0;
    
    FVL(mi) = FVL(mi+1)+Li(mi)+Wi(mi)+FIi(mi)+Reac2(mi);
    BM3(mi) = BM3(mi+1)+(FVL(mi)+Reac2(mi))*ls(mi)-Lcp(mi)+(Wi(mi)+FIi(mi))*xcgi(mi);
end
% BM2(length(masses)+2) = 0;

d = 7;
Aref = pi*(d/12/2)^2;
Lnose2(1) = 0;
Lbody2(1) = 0;
Ltail2(1) = 0;
Lfins2(1) = 0;
cgs(length(masses)+1,2) = 0;
% FVL(1) = 0;
% BM4(1) = 0;
FVL(length(masses)+1) = 0;
BM4(length(masses)+1) = 0;
% FVL(1:2) = 0;
% BM4(1:2) = 0;
% masses(1) = 0;
% distances(2) = 0;
Reac2 = 0;
r = r0;
aN2 = ((-0.5*aoa(I_maxQCNa)*sum(CNa(I_maxQCNa,:))*Aref*rho(I_maxQCNa)*norm(vrw(I_maxQCNa,:))^2+...
    Reac1(I_maxQCNa,1))/(w(I_maxQCNa)/gs(I_maxQCNa)));
for mi = length(masses):-1:1
    ks(mi) = 1/2;
    ls(mi) = distances(mi+1);
    if mi==1
        ks(mi) = 1/4;
    end
    if mi==length(masses)
        ks(mi) = 1-(1+2*(r/d)+3*(r/d)^2)/(4*(1+(r/d)+(r/d)^2));
    end
    cgs(mi,3) = ks(mi)*ls(mi)+sum(distances(mi+2:end));
    Ii(mi) = (masses(mi)*gs(I_maxQCNa)/12)*(3*(r/12)^2+(ls(mi)/12)^2);
    Iz(mi) = (masses(mi)*gs(I_maxQCNa)/2)*(r/12)^2;
end
Id = sum(Ii(1:7)+masses(1:7)*gs(I_maxQCNa).*((Rl-cg(I_maxQCNa))-transpose(cgs(1:7,3)/12)).^2);
Id2 = sum(Iz(1:7));

% dadt = Q(I_maxQCNa)*Aref*aoa(I_maxQCNa)/Id2*sum(CNa(I_maxQCNa,1:4).*...
%     (-cg(I_maxQCNa)+(CPi(I_maxQCNa,1:4))))-Reac1(I_maxQCNa,1)*((Rl-RB/12)-cg(I_maxQCNa))/Id2;
% for mi = length(masses):-1:1
%     if mi==4
%         Lbody2(mi) = CNa(I_maxQCNa,2)*aoai(I_maxQCNa,1)/2*Aref*rho(I_maxQCNa)*...
%             vrwrot(I_maxQCNa,2)^2;
%     else
%         Lbody2(mi) = 0;
%     end
%     if mi==1 %1
%         Lnose2(mi) = CNa(I_maxQCNa,1)*aoai(I_maxQCNa,1)/2*Aref*rho(I_maxQCNa)*...
%             vrwrot(I_maxQCNa,1)^2;
% %         Lnose2(mi) = 0;
%         Reac2 = Reac1(I_maxQCNa,1);
%     else
%         Lnose2(mi) = 0;
%         Reac2 = 0;
%     end
%     if mi==length(masses)
%         Ltail2(mi) = CNa(I_maxQCNa,3)*aoai(I_maxQCNa,3)/2*Aref*rho(I_maxQCNa)*vrwrot(I_maxQCNa,3)^2;
%         Lfins2(mi) = CNa(I_maxQCNa,4)*aoai(I_maxQCNa,4)/2*Aref*rho(I_maxQCNa)*vrwrot(I_maxQCNa,4)^2;
%     else
%         Ltail2(mi) = 0;
%         Lfins2(mi) = 0;
%     end
%     FVL(mi) = (Lnose2(mi)+Lbody2(mi)+Ltail2(mi)+Lfins2(mi))+FVL(mi+1)+...
%         masses(mi)*(-aN(I_maxQCNa)+((Rl-cg(I_maxQCNa))-cgs(mi,3)/12)*dadt)+...
%         masses(mi)*gs(I_maxQCNa)*sin(p(I_maxQCNa));
% %     FVL(mi+1) = (Lnose2(mi)+Lbody2(mi)+Ltail2(mi)+Lfins2(mi))+FVL(mi)+...
% %         masses(mi)*(aN2-(12*(cg(I_maxQCNa)-RB/12)-cgs(mi,3))/12*alpha(I_maxQCNa));
%     BM4(mi) = BM4(mi+1)+(cgs(mi,3)-sum(distances(mi+2:end)))*FVL(mi+1)/12+...
%         (sum(distances(mi+1:end))-cgs(mi,3))*FVL(mi)/12-...
%         Lnose2(mi)*(cgs(mi,3)-(12*(Rl-CPi(I_maxQCNa,1))))/12-...
%         Lbody2(mi)*(cgs(mi,3)-(12*(Rl-CPi(I_maxQCNa,2))))/12-...
%         Ltail2(mi)*(cgs(mi,3)-(12*(Rl-CPi(I_maxQCNa,3))))/12-...
%         Lfins2(mi)*(cgs(mi,3)-(12*(Rl-CPi(I_maxQCNa,4))))/12-...
%         Ii(mi)*dadt;
% %     BM4(mi+1) = BM4(mi)+(cgs(mi,3)-sum(distances(1:mi)))*FVL(mi)/12+...
% %         (ls(mi)+sum(distances(1:mi))-cgs(mi,3))*FVL(mi+1)/12+...
% %         Lnose2(mi)*(cgs(mi,3)-(12*(CPi(I_maxQCNa,1)-RB/12)))/12+...
% %         Lbody2(mi)*(cgs(mi,3)-(12*(CPi(I_maxQCNa,2)-RB/12)))/12+...
% %         Ltail2(mi)*(cgs(mi,3)-(12*(CPi(I_maxQCNa,3)-RB/12)))/12+...
% %         Lfins2(mi)*(cgs(mi,3)-(12*(CPi(I_maxQCNa,4)-RB/12)))/12+...
% %         Ii(mi)*dadt-Reac2*(cgs(mi,3)-distances(mi))/12;
% end
% 
Feq = abs(FAL2)+2*abs(BM3)/(r0/12);
% Feq2 = abs(FAL2)+2*abs(BM2)/(r/12);
% Feq3 = abs(FAL(2:end))+2*abs(BM2(1:end-1))/(r/12);

%% Plotting
% Plot altitude vs. time
% Plot velocity vs. time
% Plot acceleration vs. time

figure
%plot(ts(1:itf),ys(1:length(1:itf)));
%plot(ts,x(1:length(ts)));
%plot(ts(1:length(x(:,1))),x(:,1));
% plot(ts(1:itf),x(1:itf,1));
% xlabel('Time')
% ylabel('Horizontal Displacement')
% figure
%plot(ts,z(1:length(ts)));
%plot(ts(1:length(x(:,2))),x(:,2));
plot(ts(1:itf),zs(1:itf));
xlabel('Time')
ylabel('Altitude')
% figure
% %plot(ts(1:itf),vs(1:length(1:itf)));
% %plot(ts,vx(1:length(ts)));
% %plot(ts(1:length(v(:,1))),v(:,1));
% plot(ts(1:itf),v(1:itf,1));
% xlabel('Time')
% ylabel('Horizontal Velocity')
% figure
%plot(ts,vz(1:length(ts)));
%plot(ts(1:length(v(:,2))),v(:,2));
% plot(ts(1:itf),v(1:itf,2));
% xlabel('Time')
% ylabel('Vertical Velocity')
% figure
%plot(ts(1:itf-1),as(1:length(1:itf-1)));
%plot(ts(1:length(a(:,1))),a(:,1));
% plot(ts(1:itf),a(1:itf,1));
% xlabel('Time')
% ylabel('Horizontal Accceleration')
% figure
% %plot(ts(1:length(a(:,2))),a(:,2));
% plot(ts(1:itf),a(1:itf,2));
% xlabel('Time')
% ylabel('Vertical Acceleration')
figure
plot(x(1:length(zs),1),zs);
xlabel('Westing (ft)');
ylabel('Altitude (ft)');
%%axis equal
figure
plot3(x(1:length(zs),1),x(1:length(zs),2),zs);
xlabel('Westing (ft)');
ylabel('Southing (ft)');
zlabel('Altitude (ft)');
grid on;
%%axis equal
figure
%plot(ts(1:length(aoa)),aoa.*180/pi);
plot(ts(OTRi:Iapogee),aoa(OTRi:Iapogee).*180/pi);
xlabel('Time');
ylabel('aoa');
% figure
% plot(ts(1:length(v)),v);
% xlabel('Time')
% ylabel('Velocity')
% figure
% plot(ts(1:length(az)),(az.^2+ax.^2).^.5);
% xlabel('Time')
% ylabel('Acceleration')
figure
hold on
% plot(mach(1:floor(tb/dt-1)),cd8(1:floor(tb/dt-1)));
% plot(mach(floor(tb/dt-1):Iapogee-1),cd8(floor(tb/dt-1):Iapogee-1));
% plot(mach(OTRi:floor(tb/dt-1)),cd8(OTRi:floor(tb/dt-1)));
plot(mach(OTRi:floor(tb/dt-1)),cd8(OTRi:floor(tb/dt-1)));
plot(mach(OTRi:floor(tb/dt-1)),cd(OTRi:floor(tb/dt-1)));
% plot(mach(OTRi:floor(tb/dt-1)),cd(OTRi:floor(tb/dt-1))-CDfpro(OTRi:floor(tb/dt-1)));
xlabel('Mach');
ylabel('CD');
% legend('With Protuberance','Without Protuberance');
% hold on
% plot(mach(1:floor(tb/dt-1)),cd(1:floor(tb/dt-1)));
% figure
% plot(mach(800:floor(tb/dt-1)),fb(800:floor(tb/dt-1)));
% hold on
% plot(mach(floor(tb/dt-1):2800),fb(floor(tb/dt-1):2800));
% figure
%plot(ts(1:length(vwgturb(:,1))),vwgturb(:,1));
% plot(ts(1:10:5000),vwgturb(1:10:5000));
% xlabel('Time');
% ylabel('Horizontal Turbulence Velocity');
% figure
%plot(ts(1:length(vwgturb(:,1))),vwgturb(:,1));
% plot(ts(1:10:5000),vwg(1:10:5000,1));
% xlabel('Time');
% ylabel('Horizontal Wind Velocity');
% figure
%plot(ts(1:length(vwgturb(:,1))),vwgturb(:,1));
% plot(ts(1:10:5000),vwgturb(1:10:5000,2));
% xlabel('Time');
% ylabel('Vertical Turbulence Velocity');
% figure
%plot(ts(1:length(vwgturb(:,1))),vwgturb(:,1));
% plot(ts(1:10:5000),vwg(1:10:5000,2));
% xlabel('Time');
% ylabel('Vertical Wind Velocity');
% figure
% plot(zs,tempr);
% xlabel('Altitude');
% ylabel('Temperature');
% hold on
% plot(zs,temp-459.67);
% figure
% plot(zs,Pr/144);
% xlabel('Altitude');
% ylabel('Pressure');
% hold on
% plot(zs,P/144);
% figure
% plot(mach(OTRi:tb*100),CDfb(OTRi:tb*100));
% plot(mach(OTRi:tb*100),CDf(OTRi:tb*100));
% xlabel('Mach');
% ylabel('Drag');
% hold on
% plot(mach(OTRi:tb*100),CDff(OTRi:tb*100));
% plot(mach(OTRi:tb*100),CDb(OTRi:tb*100));
% plot(mach(OTRi:tb*100),CDt(OTRi:tb*100));
% plot(mach(OTRi:tb*100),CDs(OTRi:tb*100));
% plot(mach(OTRi:tb*100),CDabs(OTRi:tb*100));
% plot(mach(OTRi:tb*100),CDafs(OTRi:tb*100));
figure
plot(ts(OTRi:Iapogee),stab(OTRi:Iapogee));
xlabel('Time');
ylabel('Stability');
[MaxMach, I_MaxMach] = max(mach);
figure
plot(mach(OTRi:I_MaxMach),stab(OTRi:I_MaxMach));
xlabel('Mach');
ylabel('Stability');
% figure
% plot(ts(OTRi:Iapogee),vgust(OTRi:Iapogee));
% figure
% plot(ts(OTRi:itf),cd8(OTRi:itf));
% % hold on
% % plot(ts(OTRi:itf),cd(OTRi:itf));
% xlabel('Time');
% ylabel('CD');\
% lengths(1) = 0;
% for dDi = 1:1:length(distances)-1
%     lengths(dDi+1) = lengths(dDi)+distances(dDi+1);
% end
% figure
% plot(lengths(1:length(BM3)),12*abs(BM3));
% xlabel('Distance From Nose Tip (Inches)');
% ylabel('Bending Moment (Lbf*in)');
% figure
% plot(x(1:length(rhod),3),rhod);
% hold on
% plot(x(:,3),rho);
% xlabel('Altitude');
% ylabel('Density');
% legend('dry','wet');
figure
plot(x(:,3),rho);
xlabel('Altitude');
ylabel('Density');