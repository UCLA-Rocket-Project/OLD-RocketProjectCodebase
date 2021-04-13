function [FAL,BM,Feq,Id,Li] = loadsAnalysis(D,T,p,Reac,masses,distances,r,d,g,...
Rl,cg,CNa,aoai,rho,vrwrot,aA,aN,CPi,alpha)
% DfbT = norm(D(I_maxQCNa,:))*(CDfb(I_maxQCNa)+CDfpro(I_maxQCNa)+CDll(I_maxQCNa)+...
% CDabs(I_maxQCNa)/betas(I_maxQCNa))/cd8(I_maxQCNa); %Model protuberances separately later
DfbT = norm(D(1,:))*0.6; %Model protuberances separately later
% Dst = norm(D(I_maxQCNa,:))*(CDt(I_maxQCNa)+CDs(I_maxQCNa))/cd8(I_maxQCNa);
Dst = norm(D(1,:))*0.1;
Reac1(1,2) = dot(Reac(1,:),[-sin(p),cos(p)]);
Reac1(1,1) = dot(Reac(1,:),[-cos(p),-sin(p)]);

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
    Ii(mi) = (masses(mi)*g/12)*(3*(r/12)^2+(ls(mi)/12)^2);
end
Id = sum(Ii(1:7)+masses(1:7)*g.*((Rl-cg)-transpose(cgs(1:7,3)/12)).^2);

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
        Di(mi) = norm(D(1,:))*0.3;
    else
        Di(mi) = 0;
    end
    if mi == 1
        Reaci(mi) = Reac1(1,2);
    else
        Reaci(mi) = 0;
    end
    if mi == length(masses)
%         Ti(mi) = norm(T(I_maxQCNa-1,:)); %Max q-alpha
        Ti(mi) = norm(T(1,:));
    else
        Ti(mi) = 0;
    end
%     Dbff2 = norm(D(I_maxQCNa,:))*(CDff(I_maxQCNa)+CDb(I_maxQCNa)+...
%         CDafs(I_maxQCNa)/betas(I_maxQCNa))/cd8(I_maxQCNa);
    FIg(mi) = masses(mi)*(aA+g);
    FAL2(mi) = FAL2(mi+1)+Ti(mi)-Di(mi)-FIg(mi)+Reaci(mi);
%     FAL2(mi) = FAL2(mi+1)-Ti(mi)+Di(mi)-FIg(mi)-Reaci(mi);
%     FAL2(mi) = -sum(masses(mi:end))*abs(aA(I_maxQCNa))-Dfb22(mi);%+weight;
end
% FAL2(length(masses)+2) = 0;
FAL = FAL2;

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
Ltail(length(masses)) = -CNa(1,3)*aoai(1,3)/2*Aref*...
    rho*vrwrot(1,3)^2;
Lfins(length(masses)) = -CNa(1,4)*aoai(1,4)/2*Aref*...
    rho*vrwrot(1,4)^2;
cgs(length(masses)+1,2) = 0;
% BM2(length(masses)+1) = 0;
Reac2(length(masses)) = 0;
% aN2 = ((-0.5*aoa(I_maxQCNa)*sum(CNa(I_maxQCNa,:))*Aref*rho(I_maxQCNa)*norm(vrw(I_maxQCNa,:))^2+...
%     Reac1(I_maxQCNa,1))/(w(I_maxQCNa)/gs(I_maxQCNa)));
FVL(length(masses)+1) = 0;
BM3(length(masses)+1) = 0;

aN2 = aN;
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
        xCPi(mi,3) = Rl-CPi(1,3)-sum(distances(mi+2:end))/12;
        xCPi(mi,4) = Rl-CPi(1,4)-sum(distances(mi+2:end))/12;
        Lcp(mi) = Ltail(mi)*xCPi(mi,3)+Lfins(mi)*xCPi(mi,4);
    end
    if mi==4
        Lbody(mi) = -CNa(1,2)*aoai(1,2)/2*Aref*rho*...
            vrwrot(1,2)^2;
        xCPi(mi,2) = Rl-CPi(1,2)-sum(distances(mi+2:end))/12;
        Lcp(mi) = Lbody(mi)*xCPi(mi,2);
    else
        Lbody(mi) = 0;
    end
    FIi(mi) = masses(mi)*(-aN2+alpha*((Rl-cg)-cgs(mi,2)/12));
    if mi==2 
        Reac2(mi) = Reac1(1,1);
%         FIi(mi) = FIi(mi) -...
%             masses(mi-1)*(-aN2+alpha(I_maxQCNa)*((Rl-cg(I_maxQCNa))-cgs(mi-1,2)/12));
    else
        Reac2(mi) = 0;
    end
    if mi==1
        Lnose(mi) = -CNa(1,1)*aoai(1,1)/2*Aref*rho*...
            vrwrot(1,1)^2;
        xCPi(mi,1) = Rl-CPi(1,1)-sum(distances(mi+2:end))/12;
        Lcp(mi) = Lnose(mi)*xCPi(mi,1);
%         Lnose(mi) = 0;
%         Lcp(mi) = 0;
%         FIi(mi) = 0;
    else
        Lnose(mi) = 0;
    end
    Wi(mi) = masses(mi)*g*sin(p);
    Li(mi) = Lnose(mi)+Lbody(mi)+Ltail(mi)+Lfins(mi);
%     Li(mi) = 0;
%     Lcp(mi) = 0;
    
    FVL(mi) = FVL(mi+1)+Li(mi)+Wi(mi)+FIi(mi)+Reac2(mi);
    BM3(mi) = BM3(mi+1)+(FVL(mi)+Reac2(mi))*ls(mi)-Lcp(mi)+(Wi(mi)+FIi(mi))*xcgi(mi);
end
% BM2(length(masses)+2) = 0;
BM = BM3;

% d = 7;
% Aref = pi*(d/12/2)^2;
% Lnose2(1) = 0;
% Lbody2(1) = 0;
% Ltail2(1) = 0;
% Lfins2(1) = 0;
% cgs(length(masses)+1,2) = 0;
% % FVL(1) = 0;
% % BM4(1) = 0;
% FVL(length(masses)+1) = 0;
% BM4(length(masses)+1) = 0;
% % FVL(1:2) = 0;
% % BM4(1:2) = 0;
% % masses(1) = 0;
% % distances(2) = 0; 
% Reac2 = 0;
% aN2 = ((-0.5*aoa(I_maxQCNa)*sum(CNa(I_maxQCNa,:))*Aref*rho(I_maxQCNa)*norm(vrw(I_maxQCNa,:))^2+...
%     Reac1(I_maxQCNa,1))/(w(I_maxQCNa)/gs(I_maxQCNa)));
% for mi = length(masses):-1:1
%     ks(mi) = 1/2;
%     ls(mi) = distances(mi+1);
%     if mi==1
%         ks(mi) = 1/4;
%     end
%     if mi==length(masses)
%         ks(mi) = 1-(1+2*(r/d)+3*(r/d)^2)/(4*(1+(r/d)+(r/d)^2));
%     end
%     cgs(mi,3) = ks(mi)*ls(mi)+sum(distances(mi+2:end));
%     Ii(mi) = (masses(mi)*gs(I_maxQCNa)/12)*(3*(r/12)^2+(ls(mi)/12)^2);
%     Iz(mi) = (masses(mi)*gs(I_maxQCNa)/2)*(r/12)^2;
% end
% Id2 = sum(Ii(1:7)+masses(1:7)*gs(I_maxQCNa).*((Rl-cg(I_maxQCNa))-transpose(cgs(1:7,3)/12)).^2);
% Idz = sum(Iz(1:7));

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
Feq = abs(FAL2)+2*abs(BM3)/(r/12);
% Feq2 = abs(FAL2)+2*abs(BM2)/(r/12);
% Feq3 = abs(FAL(2:end))+2*abs(BM2(1:end-1))/(r/12);
end