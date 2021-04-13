%% Clear Contamination From Past Runs

%clc; 
clear variables; 
% clear variables Frfs 
%clearvars -except zzz
close all;
format short

%% Initial values & Initialize arrays

load('InputParams_3DOF');

%% Simulate

i = 1;
for t = ti:dt:tf
    %% There must be forces on rocket at all times
    if i>length(F)
        break
    end
    
    %% Atmosphere Model
    z = x(i,2)-h;
    zs(i) = z;
    H = x(i,2);
    Hs(i) = H;
    if i==1
        tempr0 = 59 - 0.00356*H;
        H0 = H;
        Pr0 = 473.1*exp(1)^(1.73-0.000048*H);
    end
    [temp(i),P(i),rho(i)] = atmosphere(H,Temp0,P0,tempr0,Pr0,H0,kTemp,R);
    
    %% Wind Speed Model
    
    if i~=1
        dpdxp = dpdx(i-1);
        P1 = P(i-1);
    else
        dpdxp = 0;
        P1 = P(i);
    end
    [vwg(i,1),dpdx(i)] = windspeed(vwgx,rho(i),P(i),P1,H,z,h,v(i,2),...
    ustar(i),k,z1,LMBL,hv,fv,dpdxp);
    
    %% Gravity Model
    %g = (glat-3.086*(10^(-6))*(x(i,2)/3.28))*3.28;
    g = GM/(re+H)^2;
    m = w(i)/g;
    
    %% Wind Turbulence Model
    if z >= 19 && z <= 21
        vwg20 = vwg(i,1);
    end
    
    if turb ~= 0 && z > 21
        [vwgturb(i+1,:),vwg(i,:)] = turbulence(vwg(i,:),turb,z,s(i),nufs,vwg20,vwgturb(i,:),dt,i);
    end
        
    %% Wind Gusting Model
    vgust(i) = windgust(gust,t,tgust,dtgust,dtrise);
    vwg(i,1) = vwg(i,1) + vgust(i);
    
    %% Rocket Speeds Relative to Wind
    vrw(i,:) = v(i,:)-vwg(i,:);
    vrot(i,:) = (cgi(:)-cg(i))*omega(i);
    vrotx(i,:) = vrot(i,:)*cos(p(i));
    vrotz(i,:) = vrot(i,:)*sin(p(i));
    vrwrot(i,:) = ((vrotx(i,:)+vrw(i,1)).^2+(vrotz(i,:)+vrw(i,2)).^2).^0.5;
    
    %% 0 aoa Drag Coefficient Model
    vs(i) = (gamma*R*(temp(i)+459.69))^(1/2);
    mach(i) = norm(vrw(i,:))/vs(i);
    if i==1
        CDb6 = 0;
    end
    [CD(i),CDf(i),CDb(i),CDt(i),CDs(i),CDfb(i),CDff(i),CDfpro(i),CDll(i),CDb6] = ...
    getCD0(temp(i),mu0,temp0,rho(i),vrwrot(i,:),vrw(i,:),mach(i),Lc,tr,ft,Cr,Ct,S,d,Rl,n,...
    Afe2,machD,machF,dCDmax,aPro,Apro,Spro,aLL,All,Sll,Kb,dr,nb,t,tb,CDb6);

    %% Compressibility Corrections
        
    if mach(i)>0.9 && mach(i)<1.1
        beta = sqrt(1-0.9^2);
%         CNaf = CNa(1,4)/beta;
%         CNaf = 2*pi*Afin*((Afe/2)/Aref)/(2+sqrt(4+beta*Afin/cos(thetaLE)));
    elseif mach(i) <= 0.9
        beta = sqrt(1-mach(i)^2);
%         CNaf = CNa(1,4)/beta;
%         CNaf = 2*pi*Afin*((Afe/2)/Aref)/(2+sqrt(4+beta*Afin/cos(thetaLE)));
    else
        beta = sqrt(mach(i)^2-1);
%         CNaf = CNa(1,4)/beta;
%         C3 = (gamma*mach(i)^4+(mach(i)^2-2)^2)/(2*(mach(i)^2-1)^1.5);
%         Aprime = 2*ft/3/vs(i);
%         CNaf = 4/beta*(1-(1-C3*Aprime)/2/beta/Afin);
%         CNaf = K_wb*(Sfin/Aref)*CNaf;
%         nosebodyCNa = msCNa*mach(i) + bsCNa;
%         CNan = 0.2*mach(i) + 2.03;
%         CNab = nosebodyCNa - CNan;
    end
    betas(i) = beta;
    
    %% Thrust Model
    
    if t>tb
        Tt(i) = 0;
    end
    if i==2
%     Tt(i) = 2*max(Tt);
    end
    if t<=tb 
        T(i,:) = (Tt(i) + (P(1)-P(i))*(pi*(dnoz/2/12)^2))*[-sin(p(i)+tma),cos(p(i)+tma)];
    end
    Tr(i) = norm(T(i,:))*sin(tma);
    
    %% Angles Relative to Wind
    
    if x(i,2)==h && v(i,2)==0
        phi(i) = 0;
    else
        phi(i) = asin((vwg(i,1)-v(i,1))/norm(vrw(i,:)));
        phirot(i,:) = asin((vwg(i,1)-v(i,1)-vrotx(i,:))./vrwrot(i,:));
    end
    aoa(i) = phi(i) - p(i);
    aoai(i,:) = phirot(i,:) - p(i);
    
    %% Induced Drag Coefficient Model
    
    [CDab,CDaf] = getCDi(aoa(i),aoadd,aoadn,Rl,Ln,Afp,dfa,kfb,kbf,Afe,dist(i),RL,d);
    CDabs(i) = CDab;
    CDafs(i) = CDaf;
    
    
    %% Calculate Center of Pressure
    
%     if mach(i)>1.1 && mach(i)<1.15
% %         beta = sqrt(mach(i)^2-1);
%     end
%     CNaf = CNa(1,4)*beta;
    if abs(aoai(i,2))<20*pi/180 % && mach(i)<1.2
        CNab = 4*K/pi*Ap/(d/12)^2*abs(aoai(i,2));
    else %if mach(i)<1.2
        CNab = 4*K/pi*Ap/(d/12)^2*20*pi/180;
    end
%     if mach(i)<1.2
    CNa(i,:) = [CNa(1,1)/beta,CNab/beta,CNa(1,3)/beta,CNa(1,4)/beta];
%     else
%         CNa(i,:) = [CNan,CNab,CNa(1,3)/beta,CNaf];
%     end

    if mach(i)>0.5
        CPi(i,4) = CPi(1,4)+(CPiFM2-CPiFM0)/(2-0.5)*(mach(i)-0.5);
    else
        CPi(i,4) = CPi(1,4);
    end
        
    CPi(i,1:3) = CPi(1,1:3);
    CP8(i) = sum(CNa(i,:).*CPi(i,:))/sum(CNa(i,:)); 

    if CP8(i)>Rl
        CP8(i)=Rl;
    elseif CP8(i)<0
        CP8(i)=0;
    end
    
    %% Calculate Lift Coefficient
    
    cl8(i) = sum(CNa(i,:))*aoa(i);
    
    %% Get Total Drag Coefficient
    
%     CDfd(i) = (Afe/Aref)*cl8(i)^2/(pi*Afin);
%     CDfd(i) = CDfd(i)/beta;
%     if mach(i) < 1.0
%         e = 0.9;
%     else
%         e = 0.6;
%     end
        
%     CDa(i) = cl8(i)^2/(pi*e*Afin)/beta;
%     CDa(i) = 5*sin(abs(aoa(i)));
    CDa(i) = 0;
    cd8(i) = CD(i)+CDa(i);
    cd(i) = cd8(i);
    cd8(i) = CD(i)+(CDabs(i)+CDafs(i))/beta;
    %cd8(i) = CD(i)+CDfd(i);

    if cd8(i) < 0.3
        cd8(i) = 0.3;
    elseif cd8(i) > 1.6
        cd8(i) = 1.6;
    end
    
    %% Calculate Drag Force
    
    if v(i,2) >= 0
        D(i,:) = cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*[sin(p(i)),-cos(p(i))];
    else
        phi(i) = -pi-phi(i);
        D(i,:) = cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*[sin(phi(i)),-cos(phi(i))];
        if v(i-1,2) >=0 && ejection == 0
            Iapogee2 = i;
        end
        if recovery == 0 && dist(i)>RL
            tff = t;
            itf = i;
            break
        end
        if zs(i) > 1500
            cd8(i) = 1.6;
            d = 24;
            Aref = pi*(d/12/2)^2;
            D(i,:) = cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*[sin(phi(i)),-cos(phi(i))];
        else
            if cd8(i-1) == 1.6
                Reac(i,2) = 800;
                ejection = 0;
            end
            cd8(i) = 2.2;
            d = 60;
            Aref = pi*(d/12/2)^2;
            D(i,:) = cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*[sin(phi(i)),-cos(phi(i))];
        end
    end
    
    %% Calculate Center of Mass and Stability
    
    if t>=tb
        w(i) = dw;
        Tt(i) = 0;
    end
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
    
    if v(i,2) >=0 && dist(i)>RL
        stab(i) = (CP8(i)-cg(i))/(d/12);
    end
    
    %% Calculate Lift Force
    
    L(i,:) = sign(aoa(i))*abs(cl8(i))/2*Aref*rho(i)*...
        (norm(vrw(i,:)))^2*[cos(p(i)),sin(p(i))];
    if v(i,2) < 0
        L(i,:) = [0,0];
    end
    
    %% Calculate Weight
    
    if t<=tb
        W(i,:) = [0,-w(i)];
    else
        W(i,:) = [0,-dw];
    end
    
    %% Calculate Fictitious Forces
    
    C(i,:) = 2*m*[Eomega*v(i,2)-x(i,1)*Eomega^2,-Eomega*v(i,1)-(H+h)*Eomega^2];
    %C(i,:) = [0,0];
    
    %% Calculate Rail Friction Force
    
    if s(i) > 0 && dist(i) < RL && t<tb
        Frfk(i) = muk*w(i)*sin(abs(p(i)))+2*muk/Lll*((norm(T(i,:))-w(i)*cos(p(i))-norm(D(i)))*...
            (r/12+rll)+w(i)*sin(abs(p(i)))*(aLL-cg(i)));
        Frk(i,:) = abs(Frfk(i))*[sin(p(i)),-cos(p(i))];
    elseif dist(i) <= 0 && t<tb
        Frfs(i) = mus*w(i)*sin(abs(p(i)))+2*mus/Lll*((norm(T(i,:))-w(i)*cos(p(i)))*...
            (r/12+rll)+w(i)*sin(abs(p(i)))*(aLL-cg(i)));
        Frs(i,:) = abs(Frfs(i))*[sin(p(i)),-cos(p(i))];
        Frfk(i) = muk*w(i)*sin(abs(p(i)))+2*muk/Lll*((norm(T(i,:))-w(i)*cos(p(i))-norm(D(i)))*...
            (r/12+rll)+w(i)*sin(abs(p(i)))*(aLL-cg(i)));
        Frk(i,:) = abs(Frfk(i))*[sin(p(i)),-cos(p(i))];
    end
    
    %% Calculate Reaction Forces and Moments due to Recovery Shock
    
    if Iapogee2 ~= 0 && ejection==0 && cd8(i) ~= 2.2
        %Shock = [1700/3,1700];
        Shock = [0,0];
        Reac(i,:) = Shock(1)*[-sign(v(i,1))*cos(p(i)),-sin(p(i))]; %Make perpendicular to axis
        Reac(i,:) = Reac(i,:)+Shock(2)*[-sin(p(i)),cos(p(i))]; %Axial
        Mreac = -sign(p(i))*Shock(1)*(cg(i)-47/12);
        ejection = 1;
    else
        Mreac = 0;
    end
    
    %% Calculate Total Force on Rocket
    
    if dist(i) > RL || v(i,2) < 0 % dist(i) > RL-Rl+cg(i) or dist(i) > RL for OTRS?
        F(i,:) = T(i,:)+W(i,:)+D(i,:)+L(i,:)-C(i,:);%+Reac(i,:);
        if v(i,2)>-1
            AF(i) = dot(F(i,:),[-sin(p(i)),cos(p(i))]);
            NF(i) = dot(F(i,:),[-cos(p(i)),-sin(p(i))]);
        end
        if Reac(i,2)==800
            I_main = i;
            AF(i) = norm(W(i,:))-norm(D(i,:))-norm(Reac(i,:));
            NF(i) = 0;
        end
    elseif ((norm(T(i,:))-w(i)*cos(p(i)))>abs(Frfs(end)) || s(i) > 0) && t<tb  %Moving on rail
        Reac(i,:) = [norm(W(i,:))*sin(p(i))*cos(p(i)),norm(W(i,:))*sin(p(i))^2]+Frk(i,:);
        F(i,:) = T(i,:)+W(i,:)+D(i,:)+Frk(i,:)+Reac(i,:);
        AF(i) = norm(T(i,:))-norm(D(i,:))-norm(W(i,:))*cos(p(i));
        AF2(i) = dot(F(i,:),[-sin(p(i)),cos(p(i))]);
        NF(i) = norm(W(i,:))*sin(p(i))-sign(p(i))*norm(Reac(i,:));
        OTRs = s(i);
        OTRi = i;
    elseif ((norm(T(i,:))-w(i)*cos(p(i)))<=abs(Frfs(i)) && dist(i) <= 0) && t<tb
        Reac(i,:) = -(T(i,:)+W(i,:)+D(i,:)+Frs(i,:));
        F(i,:) = T(i,:)+W(i,:)+D(i,:)+Frs(i,:)+Reac(i,:);
        L(i,:) = [0,0];
    end
    
    %% Calculate Linear Motion Variables
    
    a(i,:) = F(i,:)/m;
    v(i+1,:) = v(i,:) + a(i,:)*dt;
    x(i+1,:) = x(i,:) + v(i,:)*dt;
    dist(i+1) = norm(x(i+1,:)-x(i,:))+dist(i);
    s(i+1) = norm(v(i+1,:));
    aA(i) = AF(i)/m;
    aN(i) = NF(i)/m;
    
    %% Check for Instability
    
    if CP8(i)<cg(i) && dist(i)>RL-Rl+cg(i) && v(i,2)>0 && abs(aoa(i))<15*pi/180
%         disp(i);
        %error('Unstable');
    end
    
    %% Calculate Moment of Inertia
    
    If(i) = (1/12)*(wf(i))*(3*rtank^2 + lf(i)^2) + (wf(i))*(cg(i)-cgf(i))^2;
    Io(i) = (1/12)*(wo(i))*(3*rtank^2 + lo(i)^2) + (wo(i))*(cg(i)-cgo(i))^2;
    if t>tb
        I(i) = Id;
    else
        I(i) = Id + (dw)*(cg(i)-cgd)^2 + If(i) + Io(i);
    end
    if i ~= 1
        Idot(i) = (I(i)-I(i-1))/dt;
    end
    
    %% Calculate Moments 
    
    cgp(i) = (cgf(i)*wf(i)+cgo(i)*wo(i))/(wf(i)+wo(i));
    Mt(i) = -(wfr(i)/g)*omega(i)*((Rl-cg(i))^2-(cgp(i)-cg(i))^2); % Thrust Damping Moment
%     LM(i,:) = [cg(i)-CPi(i,1),cg(i)-CPi(i,2),CPi(i,3)-cg(i),CPi(i,4)-cg(i)];
%     MN(i) = Q(i)*Aref*sum(CNa(i,:).*aoai(i,:).*LM(i,:));
    MN(i) = (CP8(i) - cg(i))*(sign(aoa(i))*norm(L(i,:))); % Moment due to Lift Force
%     Cm(i) = -stab(i)*cl8(i);
%     MNcm(i) = Cm(i)*0.5*rho(i)*norm(vrw(i,:))^2*Aref*(d/12);
    %summ = CNa(i,1).*(CPi(1)-cg(i)).^2+CNa(i,2).*(CPi(2)-cg(i)).^2+...
    %      CNa(i,3).*(CPi(3)-cg(i)).^2+CNa(i,4).*(CPi(4)-cg(i)).^2;
%     thing = CNa(i,:).*(CPi(i,:)-cg(i)).^2;
%     summ = trace(CNa(i,:).*(CPi(i,:)-cg(i)).^2);
    summ = sum(CNa(i,:).*(CPi(i,:)-cg(i)).^2);
    %summ = CNa(i,:).*(CPi(i,:)-cg(i)).^2;
    Md(i) = 0.5*rho(i)*norm(vrw(i,:))*Aref*omega(i)*summ; % Aerodynamic Damping
%     Md(i) = 0;
    Mll(i) = 0.5*rho(i)*(norm(vrw(i,:)))^2*... % Launch Lug Moment
        (CDll(i))*Aref*(r/12); %+-Depends on which side lug is on
    Mpro(i) = 0.5*rho(i)*(norm(vrw(i,:)))^2*... % Fairing Moment
        (CDfpro(i))*Aref*(r/12); %+-Depends on protuberance side
%     Mpro(i) = 0;
    MT(i) = -(Rl-cg(i))*Tr(i); % Thrust Misalignment Moment
    M(i) = MN(i)+Mt(i)+Md(i)+Mll(i)+Mpro(i)+Mreac+MT(i); % Total Moment
%     M(i) = MN(i)+Mt(i)+Md(i);

    %% Calculate Rotational Motion Variables

    alpha(i)=(M(i)-Idot(i)*omega(i))/(I(i));
    if dist(i)<RL %|| v(i,2)<0
        alpha(i) = 0;
    end
    omega(i+1) = omega(i) + alpha(i)*dt;
    p(i+1) = p(i) + omega(i)*dt;
    
    %% Miscellaneous Values
    gs(i) = g;
    ms(i) = m;
    Q(i) = 0.5*rho(i)*(norm(vrw(i,:)))^2;
    Qa(i) = Q(i)*aoa(i);
    QCNa(i) = Q(i)*sum(abs(CNa(i,:)));

    %% End Condition
    
    if z < 0 && t > tb && v(i,2) < 0
        tff = t;
        itf = i;
        break
    end
    
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
% I_maxQCNa = Iapogee2;
% I_maxQCNa = I_main;
OTR_Tavg = sum(T(1:OTRi,2))/OTRi;
deployment_v = abs(v(Iapogee2,1));
fprintf("Apogee: %f\nDeployment Velocity: %f\nOTR Speed: %f\n",apogee,deployment_v,OTRs);

m3 = m3+wo(I_maxQCNa);
m5 = m5+wf(I_maxQCNa);
masses = [m1,m2,m3,m4,m5,m6,m7]/gs(I_maxQCNa);

[FAL,BM,Feq,Id,Li] = loadsAnalysis(D(I_maxQCNa,:),T(I_maxQCNa,:),p(I_maxQCNa),Reac(I_maxQCNa,:),...
    masses,distances,r,d,gs(I_maxQCNa),Rl,cg(I_maxQCNa),CNa(I_maxQCNa,:),aoai(I_maxQCNa,:),...
    rho(I_maxQCNa),vrwrot(I_maxQCNa,:),aA(I_maxQCNa),aN(I_maxQCNa),CPi(I_maxQCNa,:),alpha(I_maxQCNa));

%% Plotting

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
xlabel('Horizontal Displacement');
ylabel('Altitude');
axis equal
figure
%plot(ts(1:length(aoa)),aoa.*180/pi);
plot(ts(OTRi:Iapogee),abs(aoa(OTRi:Iapogee)).*180/pi);
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
% figure
% plot(ts(OTRi:Iapogee),vgust(OTRi:Iapogee));
% figure
% plot(ts(OTRi:itf),cd8(OTRi:itf));
% % hold on
% % plot(ts(OTRi:itf),cd(OTRi:itf));
% xlabel('Time');
% ylabel('CD');
lengths(1) = 0;
for dDi = 1:1:length(distances)-1
    lengths(dDi+1) = lengths(dDi)+distances(dDi+1);
end
figure
plot(lengths(1:length(BM)),12*abs(BM));
xlabel('Distance From Nose Tip (Inches)');
ylabel('Bending Moment (Lbf*in)');