function [CD,CDf,CDb,CDt,CDs,CDfb,CDff,CDfpro,CDll,CDb6] = getCD0(T,mu0,temp0,rho,vrwrot,vrw,mach,...
Lc,tr,ft,Cr,Ct,S,d,Rl,n,Afe2,machD,machF,dCDmax,aPro,Apro,Spro,aLL,All,Sll,Kb,dr,nb,t,tb,CDb6)

    temp = T+459.67;
    mu = mu0*((temp0+198.72)/(temp+198.72))*(temp/temp0)^1.5;
    Rnc(1,1) = rho*vrwrot(1,2)*Lc(1)/mu*(1+0.0283*mach-0.043*mach^2+...
        0.2107*mach^3-0.03829*mach^4+0.002709*mach^5);
    Rnc(1,2) = rho*vrwrot(1,4)*Lc(2)/mu*(1+0.0283*mach-0.043*mach^2+...
        0.2107*mach^3-0.03829*mach^4+0.002709*mach^5);
    Rnc(1,3) = rho*norm(vrw(1,:))*Lc(3)/mu*(1+0.0283*mach-0.043*mach^2+...
        0.2107*mach^3-0.03829*mach^4+0.002709*mach^5);
    Rnc(1,4) = rho*norm(vrw(1,:))*Lc(4)/mu*(1+0.0283*mach-0.043*mach^2+...
        0.2107*mach^3-0.03829*mach^4+0.002709*mach^5);
    
%     B(i,1) = Rec*(0.074./Re(i,1)^.2-1.328/sqrt(Re(i,1)));
%     B(i,2) = Rec*(0.074./Re(i,2)^.2-1.328/sqrt(Re(i,2)));
    
%     for j = 1:2
%         if Rnc(i,j)<=Rec
%             Cf(i,j) = 1.328/sqrt(Rnc(i,j));
%         else
%             Cf(i,j) = 0.074/Rnc(i,j)^.2 - B(i,j)/Rnc(i,j);
%         end
%     end
    
    Cfi(1,:) = 0.037036*Rnc(1,:).^(-0.155079);
    Cf(1,:) = Cfi(1,:)*(1+0.00798*mach-0.1813*mach^2+0.0632*mach^3-...
        0.00933*mach^4+0.000549*mach^5);
    CfiT(1,:) = (1.89+1.62*log10(Lc(:)/0.00025)).^(-2.5);
    CfT(1,:) = CfiT(1,:)/(1+0.2044*mach^2);
    if Cf(1,1)<=CfT(1,1)
        Cf(1,1) = CfT(1,1);
    end
    if Cf(1,2)<=CfT(1,2)
        Cf(1,2) = CfT(1,2);
    end
    if Cf(1,3)<=CfT(1,3)
        Cf(1,3) = CfT(1,3);
    end
    if Cf(1,4)<=CfT(1,4)
        Cf(1,4) = CfT(1,4);
    end
    SB = pi*(d/12)*Rl*0.8;
    CDfb = Cf(1,1)*(1+60/(Lc(1)/(d/12))^3+0.0025*Lc(1)/(d/12))*4*SB/pi/(d/12)^2;
    
    Rn = rho*vrwrot(1,4)*Lc(2)/mu;
    Cfg = Cf(1,2)*(log10(Rn))^2.6/(tr^2-1)*(tr^2*(log10(Rn*tr))^(-2.6)-...
        (log10(Rn))^(-2.6)+0.5646*(tr^2*(log10(Rn*tr))^(-3.6)-...
        (log10(Rn))^(-3.6)));
    CDff = Cfg*(1+60*(ft/Cr)^4+0.8*(1+5*(Cr*0.3/Cr)^2)*(ft/Cr))*4*n*Afe2/pi/(d/12)^2;
    
    X = (mach-machD)/(machF-machD);
    f = -8.3474*X^5+24.543*X^4-24.946*X^3+8.6321*X^2+1.1195*X;
    if mach<=machF && mach>=machD
        CDt = dCDmax*f;
        CDs = 0;
        dCDpro = 0.01*f;
    elseif mach > machF
        CDt = 0;
%         CDs(i) = dCDmax;
        CDs = dCDmax-0.1*abs(machF-mach);
        dCDpro = 0.01;
    else
        CDt = 0;
        CDs = 0;
        dCDpro = 0;
    end
    
    Cfpro = 0.8151*Cf(1,3)*((aPro/12)/Lc(3))^(-0.1243);
    if Apro > 0
        CDfpro = Cfpro*(1+1.798*(sqrt(Apro)/Lc(3))^1.5)*4*Spro/pi/(d/12)^2;
        CDfpro = CDfpro+dCDpro;
%         CDfpro(i) = 2*CDfpro(i);
    else
        CDfpro = 0;
    end
    %CDfpro(i) = 0;
    
    Cfll = 0.8151*Cf(1,4)*((aLL/12)/Lc(4))^(-0.1243);
    CDll = Cfll*(1+1.798*(sqrt(All)/Lc(4))^1.5)*4*Sll/pi/(d/12)^2;
    CDll = CDll+dCDpro;

    if mach<0.78
        Ke = 0.00038;
    elseif mach>=0.78 && mach<=1.04
        Ke = -0.4501*mach^4+1.5954*mach^3-2.1062*mach^2+1.2288*mach-0.26717;
    elseif mach>1.04
        Ke = 0.0002*mach^2-0.0012*mach+0.0018;
    else
        error('mach number fail');
    end
    CDfe = Ke*4*((SB+Afe2)+4*S*(Cr+Ct))/pi/(d/12)^2;
    CDfe = 0;
%     CDfe(i) = 0.5*CDfe(i);

%     if mach(i)>0.9 && mach(i)<1.1
%         beta = sqrt(1-0.9^2);
%         CNaf = CNa(1,4)/beta;
% %         CNaf = 2*pi*Afin*((Afe/2)/Aref)/(2+sqrt(4+beta*Afin/cos(thetaLE)));
%     elseif mach(i) <= 0.9
%         beta = sqrt(1-mach(i)^2);
%         CNaf = CNa(1,4)/beta;
% %         CNaf = 2*pi*Afin*((Afe/2)/Aref)/(2+sqrt(4+beta*Afin/cos(thetaLE)));
%     else
%         beta = sqrt(mach(i)^2-1);
%         CNaf = CNa(1,4)/beta;
% %         C3 = (gamma*mach(i)^4+(mach(i)^2-2)^2)/(2*(mach(i)^2-1)^1.5);
% %         Aprime = 2*ft/3/vs(i);
% %         CNaf = 4/beta*(1-(1-C3*Aprime)/2/beta/Afin);
% %         CNaf = K_wb*(Sfin/Aref)*CNaf;
% %         nosebodyCNa = msCNa*mach(i) + bsCNa;
% %         CNan = 0.2*mach(i) + 2.03;
% %         CNab = nosebodyCNa - CNan;
%     end
%     betas(i) = beta;
    CDf = CDfb+1.04*CDff+1.04*CDfpro+CDfe+1.04*CDll;
%     CDf(i) = CDf(i)/beta;
    
    if mach <= 0.6 && t<tb
       CDf6 = CDfb+1.04*CDff+1.04*CDfpro+CDfe+1.04*CDll;
       CDb6 = Kb*((dr/(d/12))^nb)/sqrt(CDf6);
    end
    
    if mach<=0.6
        CDb = Kb*((dr/(d/12))^nb)/sqrt(CDf);
    elseif mach>0.6
        if mach<=1
            fb = 1+215.8*(mach-0.6)^6;
        elseif mach<=2
            fb = 2.0881*(mach-1)^3-3.7938*(mach-1)^2+1.4618*(mach-1)+1.883917;
        elseif mach>2
            fb = 0.297*(mach-1)^3-0.7937*(mach-1)^2-0.1115*(mach-1)+1.64006;
        end
        CDb = CDb6*fb;
    end
    if t<=tb
%         CDb(i) = CDb(i)*0.8;
%         CDb(i) = 0;
    end
    
%     if t>tb
%         Tt(i) = 0;
%     end
%     if it==2
%         Tt(i) = 2*max(Tt);
%     end
%     if t<=tb 
%         T(i,:) = (Tt(i) + (P(1)-P(i))*(pi*(dnoz/2/12)^2))*[-sin(p(i)+tma),cos(p(i)+tma)];
%     end
%     Tr(i) = norm(T(i,:))*sin(tma);
%     Tx(i) = -Tt(i)*sin(p(i));
%     Tz(i) = Tt(i)*cos(p(i));

%     CDLEf(i) = n*S*ft/Aref*(((1+(gamma-1)/2*(mach(i)*cos(thetaLE))^2)^(gamma/(gamma-1))-1)...
%         /(gamma/2*(mach(i)*cos(thetaLE))^2));  %Too high
    
    CD = CDf+CDb+CDt+CDs;
end