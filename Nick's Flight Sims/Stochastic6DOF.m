clear variables;
% close all;
clc;

%% Load Map of Launch Site
img_far = imread('image1.png');
img_h_mi = 12; img_w_mi = 20;
img_h = size(img_far,1); img_w = size(img_far,2);
figure;
hold on;
RI = imref2d(size(img_far));
RI.XWorldLimits = [-img_w_mi/2 img_w_mi/2];
RI.YWorldLimits = [-img_h_mi/2 img_h_mi/2];
imshow(img_far,RI);

%% Run the Trajectory Simulations

for si = 1:1:25 % First loop is mean
clearvars -except si xs apogees OTRSs deployment_vs drift maxS
%% Initial values

recovery = 1; % 0 for no recovery / 1 for recovery
turb = 0; % 0 for none / 1 for light / 2 for medium / 3 for severe

h = 2000;
if si==1
    Temp0 = 79.5;
    kTemp = 10000; % Seems to be the most reasonable "tending to atmospheric model" rate
    P0 = 14.56*144;
    hr0 = 0.32;
else
    Temp0 = 79.5*normrnd(1,.14);
    kTemp = 10000*normrnd(1,.1); % Seems to be the most reasonable "tending to atmospheric model" rate
    P0 = 14.56*144*normrnd(1,0.005);
    hr0 = 0.32*normrnd(1,0);
end
%g = 32.174;
lat = 34.0522*pi/180;
glat = 9.7803*((1+0.001932*(sin(lat))^2)/(sqrt(1-0.006694*(sin(lat))^2)));
%g = (glat-3.086*(10^(-6))*(h/3.28))*3.28;
GM = 1.4077*10^16;
re = 20.902*10^6;
g = GM/(re+h)^2;

Eomega = 7.2921*10^(-5);

% vwgx = -5*5280/3600;
% vwgx = -27;
% vwgg = -15.2*normrnd(1,0.658/2);
% vwg = 0;
vwgd = 7*pi/15; % CCW from due west
if si==1
    vwgg = -5;
    vwgx = -vwgg*cos(vwgd)*5280/3600;
    vwgy = vwgg*sin(vwgd)*5280/3600;
else
    vwgg = -5*normrnd(1,0.658/2);
    vwgx = -vwgg*cos(vwgd)*5280/3600;
    vwgy = vwgg*sin(vwgd)*5280/3600;
end
k = 0.4;
dref = 0.000;
z1 = 0.1/100;
atm_stability = 1; % 1=neutral, 2=stable, 3=unstable

dw = 124;
load('Hot_Fire_5_1_1_Thrust_Data');
T0(:,1) = time_511;
T0(:,2) = thrust_511;
% T0 = readmatrix('Odyssey3Thrust.xlsx');
% T0 = readmatrix('PropSimThrust.xlsx');
Imp = trapz(T0(:,1),T0(:,2));
Impfactor = 9200/Imp-0.037;
% T0 = T0*Impfactor;
% Imp = trapz(T0(:,1),T0(:,2));
% Isp = 200;
wp = 50;
OF = 1.5;
wfi = wp/(OF+1);
woi = wfi*OF;
ww = dw+wfi+woi;
Isp = Imp/(ww-dw);
cgfi = 124.43/12*1.0285;
cgoi = 96.18/12*1.0285;
lfi = 22.25/12;
loi = 25.07/12;
xfuel0 = cgfi-lfi/2;
xox0 = cgoi-loi/2;
dtank = 7/12;
rtank = dtank/2;
% dnoz = 4;
dnoz = 2.81;
cgd = 111/12;
if si==1
    Ecg = 1;
    RLa = -5*pi/180;
    RL = 53.5;
else
    Ecg = normrnd(1,0.005);
    RLa = -5*pi/180*normrnd(1,0.2);
    RL = 53.5*normrnd(1,0.03);
end
p0 = RLa;
RLv = [-RL*sin(RLa),RL*cos(RLa)];
df = (wfi/g)/(lfi*pi*rtank^2);
dox = (woi/g)/(loi*pi*rtank^2);
ejection = 0;

ti = 0;
tb = T0(end,1);
% tb = round(find(T0(:,2),1,'last')/1000,2); %PropSim
tf = 750;
dt = 0.01;
t_steps = tf/dt;

Tt = interp1(T0(:,1),T0(:,2),ti:dt:tb);
% if dt==0.001
%     Tt = Tt(121:end); %If Odyssey Motor is used
%     tb = tb-0.12;
% elseif dt==0.01
%     Tt = Tt(13:end); %If Odyssey Motor is used
%     tb = tb-0.12;
% end
% Tx = zeros(t_steps,1);
% Tz = zeros(t_steps,1);
T = zeros(t_steps,3);

mach_cd = readmatrix('varCD.CSV');
ma = mach_cd(:,1);
cdv = mach_cd(:,6);
vs = zeros(t_steps,1);
gamma = 1.4;
R = 1717;
mach = zeros(t_steps,1);
cd8 = zeros(t_steps,1);
cl8 = zeros(t_steps,1);

d = 7;
Aref = pi*(d/12/2)^2;

n = 4;
CNan = 2;
dr = 5.5/12;
DR = dr/2;
CNat = 2*((dr/(d/12))^2-1);
r0 = d/2;
S = 7/12;
% S = 1/12;
Cr = 14/12;
Ct = 3/12;
Lf = sqrt(S^2+(Cr/2-Ct/2)^2);
dfa = dr;
Kf = 1+(dr/2)/(S+(dr/2)); 
CNaf = Kf*((4*n*(S/(d/12))^2)/(1+sqrt(1+(2*Lf/(Cr+Ct))^2)));
Ln = 35/12*1.0285;
xn = 0.466*Ln;
LB = 125.23/12*1.0285;
xp = Ln+LB;
Lt = 15.75/12*1.0285;
xt = xp+Lt/3*((1+(1-(d/12)/dr)/(1-((d/12)/dr)^2)));
Rl = Ln+LB+Lt;
xB = Rl-Cr;
xr = Cr-Ct;
xf = xB+xr/3*(Cr+2*Ct)/(Cr+Ct)+(1/6*((Cr+Ct)-(Cr*Ct)/(Cr+Ct)));
xBody = Ln+0.5*LB;
K = 1.5;
ApB = (d/12)*LB;
Apb = Lt/2*((d/12)+dr);
Apn = (d/12)*Ln/2*1.3; % Factor of 1.3 to account for cone to ogive difference
Ap = ApB+Apb+Apn;
CNab = 0;
Tf = 0.6/12;
thetaLE = atan(xr/S);
thetaTE = atan((Cr-xr-Ct)/S);
ft = 0.6/12;
Afe2 = 2*(Cr*S-0.5*S^2*(tan(thetaLE)+tan(thetaTE)))+ft*S*0.3;
Afe = S*(Cr+Ct);
Afp = Afe+0.5*(1/12)*Cr;
Rs = (2*S+d/12)/dfa;
kfb = 0.8065*Rs^2+1.1553*Rs;
kbf = 0.1935*Rs^2+0.8174*Rs+1;
tr = Ct/Cr;
if si==1
    Ecp = 1;
    Ecl = 1;
    Ecd = 1;
else
    Ecp = normrnd(1,0.02);
    Ecl = normrnd(1,0.01);
    Ecd = normrnd(1,0.03);
end

Apf = Afe;
Apf = 0;
Apt = Ap + Apf;
CPcut = (Apn*xn+ApB*xB+Apb*xt+Apf*xf)/Apt;

Iapogee2 = 0;

Aprime = 2*ft/3/Cr;
% Afin = (2*S)^2/(Afe/2);
Sfin = (Cr+Ct)*S;
Afin = (2*S)^2/Sfin;
K_wb = 1 + (((d/12)/2)/S)^1.15;
beta2 = sqrt(3);
MAC = (Cr+Ct)/2;
CPiFM2r = ((Afin*beta2-0.67)/(2*Afin*beta2-1));
CPiFM0 = xf-xB;
CPiFM2 = (xf-xB)*CPiFM2r/0.25;

theta = rad2deg(atan(((d/12)/2)/Ln));
msCNa = 0.05525641*theta - 0.067948718;
bsCNa = -0.044717949*theta + 2.03025641;

% machD = 0.8;
machD = -0.0156*(Ln/(d/12))^2+0.136*(Ln/(d/12))+0.6817;
if (Ln/Rl)<0.2
    A = 2.4;
    b = -1.05;
else
    A = -321.94*(Ln/Rl)^2+264.07*(Ln/Rl)-36.348;
    b = 19.634*(Ln/Rl)^2-18.369*(Ln/Rl)+1.7434;
end
machF = A*(Rl/(d/12))^b+1.0275;
c = 50.676*(Ln/Rl)^2-51.734*(Ln/Rl)+15.642;
G = -2.2538*(Ln/Rl)^2+1.3108*(Ln/Rl)-1.7344;
if (Rl/(d/12))>=6
    dCDmax = c*(Rl/(d/12))^G;
else
    dCDmax = c*(6)^G;
end
% dCDmax = 0.086;
% dCDmax = 0.115;
dCDmax = 0.1956;

Rec = 5*10^5;
mu0 = 3.62*10^(-7);
temp0 = 518.7;

Kb = 0.0274*atan((Rl-Ln)/(d/12)+0.0116);
nb = 3.6542*((Rl-Ln)/(d/12))^(-0.2733);

% Apro = 6/144;
% Apro = 0.84/144;
% Apro = 0;
rpro = 0.97/12;
Apro = pi/2*rpro^2;
Lp = 68.4/12;
aPro = 75.08;
% rpro = sqrt(Apro/pi);
Spro = 2*Apro+2*pi*rpro*Lp;
% Spro = 2*2*Apro+2*2*pi*rpro/2*Lp;
% Spro = 113.03/144;

All = 1/144; %in
rll = sqrt(All/pi);
Lll = 1.72/12;
Sll = 2*All+2*pi*rll*Lll;
aLL = 80.97;
aLL2 = 150.37;
muk = 0.25; %Al on Al lubricated
mus = 0.3; %Al on Al lubricated
% muk = 0.0; %Al on Al lubricated
% mus = 0.0; %Al on Al lubricated
if p0<0
%     mus = -mus;
end

Lc = [Rl,Cr,Lp,Lll];

%Find cg formulas later
cgn = 19.95/12*1.0285;
cgb = Ln+LB/2*1.0285;
cgt = 157.63/12*1.0285;
cgfin = 164.067/12*1.0285;
cgi = [cgn, cgb, cgt, cgfin];

vwg = zeros(t_steps,3); % Wind model in loop
ustar = zeros(t_steps,1);
% ustar(:) = abs(vwgx)/3.28*k/log((RLv(2)/3.28-dref/3.28)/(z1/3.28));
cv = 0.25; % 0.1 to 0.4
fv = 2*Eomega*sin(lat);
% ustar(:) = vwgx*k/log(RLv(2)/z1);
if atm_stability == 1
    Lv = 500*3.281;
    ustar(:) = abs(vwgx)*k/log(RLv(2)/z1);
    hv = cv*ustar(1)/abs(fv);
    vA = 4.9;
    vB = 1.9;
    LMBL = hv/2*(((log(ustar(1)/abs(fv)/z1)-vB)^2+vA^2)^0.5-log(hv/z1))^(-1);
elseif atm_stability == 2 % How to find hv without ustar
    Lv = 125*3.281;
%     ustar(:) = abs(vwgx)*k/(log(RLv(2)/z1)+4.7*RLv(2)/Lv*(1-RLv(2)/2/hv));
%     hv = 0.5*sqrt(ustar(1)*Lv/abs(fv)); 
    vA = log(hv/Lv)-2.2*hv/Lv+2.9;
    vB = 3.5*hv/Lv;
    LMBL = hv/2*(((log(ustar(1)/abs(fv)/z1)-vB)^2+vA^2)^0.5-log(hv/z1)-4.7*hv/2/Lv)^(-1);
elseif atm_stability == 3
    Lv = -150*3.281;
    vx = (1-12*RLv(2)/Lv)^(1/3);
    fzl = 1.5*log((1+vx+vx^2)/3)-sqrt(3)*atan((1+2*vx)/sqrt(3))+pi/sqrt(3);
    ustar(:) = abs(vwgx)*k/(log(RLv(2)/z1)-fzl);
    hv = cv*ustar(1)/abs(fv);
    vxh = (1-12*hv/Lv)^(1/3);
    fhl = 1.5*log((1+vxh+vxh^2)/3)-sqrt(3)*atan((1+2*vxh)/sqrt(3))+pi/sqrt(3);
    vA = log(hv/abs(Lv))+log(abs(fv)*hv/ustar(1))+1.5;
    vB = k*ustar(1)/fv/hv+1.8*fv*hv/ustar(1)*exp(1)^(0.2*hv/Lv);
    LMBL = hv/2*(((log(ustar(1)/abs(fv)/z1)-vB)^2+vA^2)^0.5-log(hv/z1)+fhl)^(-1);
else 
    error('Fix atm stability');
end
% if vwgx<0
%     ustar(:) = -ustar(:);
% end

% tgust = 5.82;
tgust = 15.93;
dtgust = 3;
dtrise = 0*dtgust;
if (rand<0.1)
    gust = -30;
else
    gust = 0;
end

% tma = [0.16/180*pi,0];
if si==1
    tmap = 0;
    tma = 0*[cos(tmap),sin(tmap)];
else
    tmap = 2*pi*rand(1);
    tma = normrnd(0,0.01)*[sin(tmap),cos(tmap)];
end
% tmap = 0.7;
% tma = 0.32/180*pi*[cos(tmap),sin(tmap)];

% cant = 1/180*pi;
if si==1
    cant = 0;
else
    cant = normrnd(0,0);
end
% rf = r0/12;
% CR = 0.04;

%% Initialize arrays

ts = ti:dt:tf;

rho = zeros(t_steps,1);
%temp = zeros(t_steps,1);
%P = zeros(t_steps,1);

wfr = zeros(t_steps,1);
w = zeros(t_steps,1);
w(1) = ww;
W = zeros(t_steps,3);
wf = zeros(t_steps,1);
wf(1) = wfi;
wo = zeros(t_steps,1);
wo(1) = woi;
cgf = zeros(t_steps,1);
cgf(1) = cgfi;
cgo = zeros(t_steps,1);
cgo(1) = cgoi;
lf = zeros(t_steps,1);
lf(1) = lfi;
lo = zeros(t_steps,1);
lo(1) = loi;
cg = zeros(t_steps,1);
wfrf = zeros(t_steps,1);
wfro = zeros(t_steps,1);
cgp = zeros(t_steps,1);
Mt = zeros(t_steps,3);
stab = zeros(t_steps,1);

% Reb = zeros(t_steps,1);
% Ref = zeros(t_steps,1);
Rnc = zeros(t_steps,4);

mu = zeros(t_steps,1);
% Bb = zeros(t_steps,1);
% Bf = zeros(t_steps,1);
B = zeros(t_steps,2);
% Cfb = zeros(t_steps,1);
% Cff = zeros(t_steps,1);
Cf = zeros(t_steps,4);
% CDfb = zeros(t_steps,1);
% CDb = zeros(t_steps,1);
% CDf = zeros(t_steps,1);
% CDi = zeros(t_steps,1);
% CD0 = zeros(t_steps,4);
Cfi = zeros(t_steps,4);
CfiT = zeros(t_steps,4);
CfT = zeros(t_steps,4);
CDf6 = zeros(t_steps,1);
% CDb6 = zeros(t_steps,1);

p = zeros(t_steps,1);
p(1) = p0;
y = zeros(t_steps,1);
y(1) = 0;
r = zeros(t_steps,1);
r(1) = 0;
omega = zeros(t_steps,3);
ihat = [-sin(p(1)),0,cos(p(1))];
jhat = [0,1,0];
khat = [cos(p(1)),0,sin(p(1))];
yA = zeros(3,t_steps);
pA = zeros(3,t_steps);
rA = zeros(3,t_steps);
q = zeros(4,t_steps);
qdot = zeros(4,t_steps);
omega(1,:) = [0,0,0]; %1=y,2=p,3=r
q(:,1) = [cos(r(1)/2)*cos(p(1)/2)*cos(y(1))+sin(r(1)/2)*sin(p(1)/2)*sin(y(1));...
     sin(r(1)/2)*cos(p(1)/2)*cos(y(1))-cos(r(1)/2)*sin(p(1)/2)*sin(y(1));...
     cos(r(1)/2)*sin(p(1)/2)*cos(y(1))+sin(r(1)/2)*cos(p(1)/2)*sin(y(1));...
     cos(r(1)/2)*cos(p(1)/2)*sin(y(1))-sin(r(1)/2)*sin(p(1)/2)*cos(y(1))];
qdot(:,1) = 0.5*([0,-omega(1,1),-omega(1,2),-omega(1,3);omega(1,1),0,omega(1,3),-omega(1,2);...
    omega(1,2),-omega(1,3),0,omega(1,1);omega(1,3),omega(1,2),-omega(1,1),0]*q(:,1));
% roll0 = atan2(2*(q(1)*q(2)+q(3)*q(4)),1-2*(q(2)^2+q(3)^2));
% pitch0 = asin(2*(q(1)*q(3)-q(4)*q(2)));
% yaw0 = atan2(2*(q(1)*q(4)+q(2)*q(3)),1-2*(q(3)^2+q(4)^2));
RM = [1-2*q(3)^2-2*q(4)^2, 2*q(2)*q(3)-2*q(1)*q(4), 2*q(2)*q(4)+2*q(1)*q(3);...
      2*q(2)*q(3)+2*q(1)*q(4), 1-2*q(2)^2-2*q(4)^2, 2*q(3)*q(4)-2*q(1)*q(2);...
      2*q(2)*q(4)-2*q(1)*q(3), 2*q(3)*q(4)+2*q(1)*q(2), 1-2*q(2)^2-2*q(3)^2];
RM(1:2,:) = -RM(1:2,:);
RMrec{1} = RM;
yA(:,1) = RM*[1;0;0];
pA(:,1) = RM*[0;1;0];
rA(:,1) = RM*[0;0;1];
% rA(1,1) = -rA(1,1);
tA(:,1) = [cos(tma(2)),-sin(tma(2)),0;sin(tma(2)),cos(tma(2)),0;0,0,1]*...
    [cos(tma(1)),0,sin(tma(1));0,1,0;-sin(tma(1)),0,cos(tma(1))]*rA(:,1);

L = zeros(t_steps,3);
%Lx = zeros(t_steps,1);
%Lz = zeros(t_steps,1);
CNa = zeros(t_steps,4);
CNa(1,:) = [CNan,CNab,CNat,CNaf];
CP8 = zeros(t_steps,1);
CPi = zeros(t_steps,4);
CPi(1,:) = [xn,xBody,xt,xf];
D = zeros(t_steps,3);
MN = zeros(t_steps,3);
% Dx = zeros(t_steps,1);
% Dz = zeros(t_steps,1);
phi = zeros(t_steps,1);
F = zeros(t_steps,3);
% Fx = zeros(t_steps,1);
% Fz = zeros(t_steps,1);
C = zeros(t_steps,3);
Md = zeros(t_steps,1);
M = zeros(t_steps,3);
Frk = zeros(t_steps,3);
Frfk = zeros(t_steps,1);
Reac = zeros(t_steps,3);

Lw = zeros(t_steps,1);
Lu = zeros(t_steps,1);
xTurb = zeros(t_steps,1);
zTurb = zeros(t_steps,2);

a = zeros(t_steps,3);
% ax = zeros(t_steps,1);
% az = zeros(t_steps,1);

% vx = zeros(t_steps,1);
% vz = zeros(t_steps,1);
v = zeros(t_steps,3);
v(1,1) = 0;
v(1,2) = 0;
v(1,3) = 0;


x = zeros(t_steps,3);
%z = zeros(t_steps,1);
x(1,1) = 282;
x(1,2) = 30.5;
x(1,3) = h;

dist = zeros(t_steps,1);
dist(1) = 0;

s = zeros(t_steps,1);
s(1) = 0;

%zs = zeros(t_steps,1);

Ixx = zeros(t_steps,1);
Iyy = zeros(t_steps,1);
Izz = zeros(t_steps,1);
I = zeros(t_steps,3);
% Id = 1788.5;
Id = 1430;
Idz = 6.0012;
If = zeros(t_steps,1);
Io = zeros(t_steps,1);
alpha = zeros(t_steps,3);
Idot = zeros(t_steps,1);

aoa = zeros(t_steps,1);
aoan = zeros(t_steps,1);
aoai = zeros(t_steps,4);
vrw = zeros(t_steps,3);
% vn = zeros(t_steps,1);
% vnx = zeros(t_steps,1);
% vnz = zeros(t_steps,1);
% vrwn = zeros(t_steps,2);
vrot = zeros(t_steps,3);
vrotx = zeros(t_steps,4);
vrotz = zeros(t_steps,4);
vrwrot = zeros(t_steps,4);
phirot = zeros(t_steps,4);

vwgturb = zeros(t_steps,2);
ts_turb = dt*(0:t_steps-1);
filt = fir1(64,0.5);
nu = randn(2,t_steps); % Need randn(t_steps) for one of the models
nuf = filter(filt,1,nu);
variance = 1.0;
for turbk = 1:1:2
    scale = sqrt(variance)/std(nuf(turbk,:));
    nufs(turbk,:) = scale*nuf(turbk,:);
end

vgust = zeros(t_steps,1);

j = zeros(t_steps,1);
iters = zeros(t_steps,1);

aoad = 4:2:20;
aoadd = [0.78,0.82,0.92,0.94,0.96,0.97,0.975,0.98,0.985];
aoadn = [0.6,0.63,0.66,0.68,0.72,0.73,0.74,0.75,0.76];

AF = zeros(t_steps,1);
NF = zeros(t_steps,1);


%% Simulate

i = 1;

for t = ti:dt:tf
    if i>length(F)
%         disp(i)
        break
    end
    if t>tb
%         cdv = mach_cd(:,7);
    end
    
    H = x(i,3);
    Hs(i) = H;
    z = x(i,3)-h;
    zs(i) = z;
    
    if z < 0 && t > tb && v(i,3) < 0
        tff = t;
        itf = i;
    %     disp(i+2)
        break
    end
    
    if H < 36152 && H >= 0
        temp(i) = 59 - 0.00356*H;
        tempr(i) = temp(i);
        temp(i) = tempr(i)+(Temp0-tempr(1))*exp((Hs(1)-H)/kTemp);
        P(i) = 2116*((temp(i)+459.67)/518.6)^5.256;
        Pr(i) = P(i);
        P(i) = Pr(i)+(P0-Pr(1))*exp((Hs(1)-H)/kTemp);
    elseif H >= 36152 && H < 82854
        temp(i) = -70;
        tempr(i) = temp(i);
        temp(i) = tempr(i)+(Temp0-tempr(1))*exp((Hs(1)-H)/kTemp);
        P(i) = 473.1*exp(1)^(1.73-0.000048*H);
        Pr(i) = P(i);
        P(i) = Pr(i)+(P0-Pr(1))*exp((Hs(1)-H)/kTemp);
    else
        temp(i) = -205.05+0.00164*H;
        tempr(i) = temp(i);
        temp(i) = tempr(i)+(Temp0-tempr(1))*exp((Hs(1)-H)/kTemp);
        P(i) = 51.97*((temp(i)+459.67)/389.98)^-11.388;
        Pr(i) = P(i);
        P(i) = Pr(i)+(P0-Pr(1))*exp((Hs(1)-H)/kTemp);
    end
    hr(i) = hr0;
    rhod(i) = P(i)/(R*(temp(i)+459.67));
    rho(i) = rhod(i)*(1+hr(i))/(1+1.609*hr(i));
    
%     if dist(i)>RL-Rl+cg(i)
%         if H<6560
%             vwg(i,1) = (ustar(i)/k*log((H/3.28-dref/3.28)/(z1/3.28)))*3.28;
%             %ustar(i+1) = k*(vwg(i,1)-vwg(i-1,1))/log(((x(i,2)-h)-dref)/((x(i-1,2)-h)-dref));
%             %vwg(i,1) = vwgx;
%         elseif H>=6560
%             vwg(i,1) = vwg(i-1,1)*(H/(x(i-1,2)-h))^(1/7);
%             %break
%         end
%     end
    if x(i,3) == h
        vwg(i,1) = vwgx;
    elseif z<=hv % Implement geostrophic patterns at z > hv
        if atm_stability == 1
            vwg(i,1) = ustar(i)/k*(log(z/z1)+z/LMBL-z/hv*(z/2/LMBL));
        elseif atm_stability == 2
            vwg(i,1) = ustar(i)/k*(log(z/z1)+4.7*z/Lv*(1-z/2/hv)+z/LMBL-z/hv*(z/2/LMBL));
        elseif atm_stability == 3
            vx = (1-12*z/Lv)^(1/3);
            fzl = 1.5*log((1+vx+vx^2)/3)-sqrt(3)*atan((1+2*vx)/sqrt(3))+pi/sqrt(3);
            vwg(i,1) = ustar(i)/k*(log(z/z1)-fzl+z/LMBL-z/hv*(z/2/LMBL));
        end
        if v(i,3)>=0
            dpdx(i) = vwg(i,1)*rho(i)*fv;
            dpdxi = i;
        end
    elseif z>hv
        if vwgx == 0
            dpdx(i) = 0;
        else
            dpdx(i) = ((P(i)/P(i-1)+1)/2)*dpdx(i-1);
%             dpdx(i) = dpdx(i-1);
%             dpdx(i) = d(P(i)/P(i-1))*dpdx(i-1);
        end
        vwg(i,1) = dpdx(i)/fv/rho(i);
    end
    
    %g = (glat-3.086*(10^(-6))*(x(i,2)/3.28))*3.28;
    g = GM/(re+H)^2;
    
    m = w(i)/g;
    
    if z >= 9 && z <= 11 % or 12 or 13
        vwg10 = vwg(i,1);
    end
    if z >= 19 && z <= 21 % or 22 or 23
        vwg20 = vwg(i,1);
    end
    
    if turb ~= 0 && z > 21
        if z>=21 && z<=1000
            Lw(i) = z;
            Lu(i) = z/(0.177+0.000823*z)^1.2;
            sigw(i) = 0.1*vwg20;
            sigu(i) = sigw(i)/(0.177+0.000823*z)^0.4;
        elseif z>=2000
            Lw(i) = 1750;
            Lu(i) = 1750;
            if turb == 1  % Change to the high altitude chart later
                sigw(i) = 5;
            elseif turb == 2
                sigw(i) = 10;
            elseif turb == 3
                sigw(i) = 15;
            end
            sigu(i) = sigw(i);
        elseif z>1000 && z<2000
            Lw(i) = 1000+(z-1000)*750/1000;
            Lu(i) = Lw(i);
            if turb == 1  % Change to the high altitude chart later
                sigw(i) = 0.1*vwg20+(z-1000)*(5-0.1*vwg20)/1000;
            elseif turb == 2
                sigw(i) = 0.1*vwg20+(z-1000)*(10-0.1*vwg20)/1000;
            elseif turb == 3
                sigw(i) = 0.1*vwg20+(z-1000)*(15-0.1*vwg20)/1000;
            end
            sigu(i) = sigw(i);
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
        vwgturb(i+1,1) = (1-s(i)*(dt/100)/Lu(i))*vwgturb(i,1)+...
            sqrt(2*s(i)*(dt/100)/Lu(i))*sigu(i)/std(nufs(1,:))*nufs(1,i);
        vwgturb(i+1,2) = (1-s(i)*(dt/100)/Lu(i))*vwgturb(i,2)+...
            sqrt(2*s(i)*(dt/100)/Lu(i))*sigw(i)/std(nufs(2,:))*nufs(2,i);
        vwg(i,:) = vwg(i,:) + vwgturb(i,:);
    end
    
    if t >= tgust && t <= tgust+dtgust+2*dtrise
        if t <= tgust+dtrise
            if dtrise == 0
                vgust(i) = gust;
            else
                vgust(i) = gust/2*(1-cos(pi*(t-tgust)/dtrise));
            end
        elseif t <= tgust+dtrise+dtgust
            vgust(i) = gust;
        else
            if dtrise == 0
                vgust(i) = gust;
            else
                vgust(i) = gust/2*(1-cos(pi*(t-tgust-dtgust)/dtrise));
            end
        end    
    else
        vgust(i) = 0;
    end
    
    vwg(i,1) = vwg(i,1) + vgust(i);
    
%     vwg(i,1) = vwgx;
    vwg(i,2) = vwgy;
    
    %vrw(i) = ((vx(i)-vwg)^2+vz(i)^2)^(1/2);
    vrw(i,:) = v(i,:)-vwg(i,:);
    aoacm(i) = acos(dot((vrw(i,:)/norm(vrw(i,:))),rA(:,i)));
    if v(i,3)<0
        aoacm(i) = 0;
    end
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
    CP8(i) = sum(CNa(i,:).*CPi(i,:))/sum(CNa(i,:))*Ecp;
%     CP(i) = CP8(i);
%     CPb(i) = sum(CNa(i,1:3).*CPi(i,1:3))/sum(CNa(i,1:3));
%     CP8(i) = CP8(i)+sin(abs(aoa(i)))*(CPcut-CP(i));

    %CP8(i) = CP8(i)/beta;
    if CP8(i)>Rl
        CP8(i)=Rl;
    elseif CP8(i)<0
        CP8(i)=0;
    end
    
    if t>tb
        Tt(i) = 0;
    end
    if i==2
%         Tt(i) = 2*max(Tt);
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
        T(i,:) = (Tt(i) + (P(1)-P(i))*(pi*(dnoz/2/12)^2))*transpose(tA(:,i));
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
        cg(i) = (cgd*dw + cgf(i)*wf(i) + cgo(i)*wo(i))/(w(i))*Ecg;
    else
        cg(i) = cgd*Ecg;
        cg(i+1) = cgd*Ecg; %To cover blind spots
    end
    if v(i,3) >=0 && dist(i)>RL
        stab(i) = (CP8(i)-cg(i))/(d/12);
    end
    
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
        vrotb(i,:) = vrw(i,:)+vrot(i,:)/abs(CP8(i)-cg(i))*(cgi(2)-cg(i));
        vrotf(i,:) = vrw(i,:)+vrot(i,:)/abs(CP8(i)-cg(i))*(cgi(4)-cg(i));
%     end
    
%     vrot(i,:) = (cgi(:)-cg(i))*omega(i);
%     vrotx(i,:) = vrot(i,:)*cos(p(i));
%     vrotz(i,:) = vrot(i,:)*sin(p(i));
%     vrwrot(i,:) = ((vrotx(i,:)+vrw(i,1)).^2+(vrotz(i,:)+vrw(i,2)).^2).^0.5;
    
    mach(i) = norm(vrw(i,:))/vs(i);
    j(i) = mach(i)*100;
%     if mach(i) == 0 
%         cd(i) = cdv(1);
%     else
%         cd(i) = cdv(floor(j(i)));
%     end
    temp(i) = temp(i)+459.67;
    mu(i) = mu0*((temp0+198.72)/(temp(i)+198.72))*(temp(i)/temp0)^1.5;
    Rnc(i,1) = rho(i)*norm(vrotb(i,:))*Lc(1)/mu(i)*(1+0.0283*mach(i)-0.043*mach(i)^2+...
        0.2107*mach(i)^3-0.03829*mach(i)^4+0.002709*mach(i)^5);
    Rnc(i,2) = rho(i)*norm(vrotf(i,:))*Lc(2)/mu(i)*(1+0.0283*mach(i)-0.043*mach(i)^2+...
        0.2107*mach(i)^3-0.03829*mach(i)^4+0.002709*mach(i)^5);
    Rnc(i,3) = rho(i)*norm(vrw(i,:))*Lc(3)/mu(i)*(1+0.0283*mach(i)-0.043*mach(i)^2+...
        0.2107*mach(i)^3-0.03829*mach(i)^4+0.002709*mach(i)^5);
    Rnc(i,4) = rho(i)*norm(vrw(i,:))*Lc(4)/mu(i)*(1+0.0283*mach(i)-0.043*mach(i)^2+...
        0.2107*mach(i)^3-0.03829*mach(i)^4+0.002709*mach(i)^5);
    
%     B(i,1) = Rec*(0.074./Re(i,1)^.2-1.328/sqrt(Re(i,1)));
%     B(i,2) = Rec*(0.074./Re(i,2)^.2-1.328/sqrt(Re(i,2)));
    
%     for j = 1:2
%         if Rnc(i,j)<=Rec
%             Cf(i,j) = 1.328/sqrt(Rnc(i,j));
%         else
%             Cf(i,j) = 0.074/Rnc(i,j)^.2 - B(i,j)/Rnc(i,j);
%         end
%     end
    
    Cfi(i,:) = 0.037036*Rnc(i,:).^(-0.155079);
    Cf(i,:) = Cfi(i,:)*(1+0.00798*mach(i)-0.1813*mach(i)^2+0.0632*mach(i)^3-...
        0.00933*mach(i)^4+0.000549*mach(i)^5);
    CfiT(i,:) = (1.89+1.62*log10(Lc(:)/0.00025)).^(-2.5);
    CfT(i,:) = CfiT(i,:)/(1+0.2044*mach(i)^2);
    if Cf(i,1)<=CfT(i,1)
        Cf(i,1) = CfT(i,1);
    end
    if Cf(i,2)<=CfT(i,2)
        Cf(i,2) = CfT(i,2);
    end
    if Cf(i,3)<=CfT(i,3)
        Cf(i,3) = CfT(i,3);
    end
    if Cf(i,4)<=CfT(i,4)
        Cf(i,4) = CfT(i,4);
    end
    SB = pi*(d/12)*Rl*0.8;
    CDfb(i) = Cf(i,1)*(1+60/(Lc(1)/(d/12))^3+0.0025*Lc(1)/(d/12))*4*SB/pi/(d/12)^2;
    
    Rn(i) = rho(i)*norm(vrotf(i,:))*Lc(2)/mu(i);
    Cfg(i) = Cf(i,2)*(log10(Rn(i)))^2.6/(tr^2-1)*(tr^2*(log10(Rn(i)*tr))^(-2.6)-...
        (log10(Rn(i)))^(-2.6)+0.5646*(tr^2*(log10(Rn(i)*tr))^(-3.6)-...
        (log10(Rn(i)))^(-3.6)));
    CDff(i) = Cfg(i)*(1+60*(ft/Cr)^4+0.8*(1+5*(Cr*0.3/Cr)^2)*(ft/Cr))*4*n*Afe2/pi/(d/12)^2;
    
    X = (mach(i)-machD)/(machF-machD);
    f(i) = -8.3474*X^5+24.543*X^4-24.946*X^3+8.6321*X^2+1.1195*X;
    if mach(i)<=machF && mach(i)>=machD
        CDt(i) = dCDmax*f(i);
        CDs(i) = 0;
        dCDpro(i) = 0.01*f(i);
    elseif mach(i) > machF
        CDt(i) = 0;
%         CDs(i) = dCDmax;
        CDs(i) = dCDmax-0.1*abs(machF-mach(i));
        dCDpro(i) = 0.01;
    else
        CDt(i) = 0;
        CDs(i) = 0;
        dCDpro(i) = 0;
    end
    
    Cfpro(i) = 0.8151*Cf(i,3)*((aPro/12)/Lc(3))^(-0.1243);
    if Apro > 0
        CDfpro(i) = Cfpro(i)*(1+1.798*(sqrt(Apro)/Lc(3))^1.5)*4*Spro/pi/(d/12)^2;
        CDfpro(i) = CDfpro(i)+dCDpro(i);
    else
        CDfpro(i) = 0;
    end
    %CDfpro(i) = 0;
    
    Cfll(i) = 0.8151*Cf(i,4)*((aLL/12)/Lll)^(-0.1243);
    CDll(i) = Cfll(i)*(1+1.798*(sqrt(All)/Lc(4))^1.5)*4*Sll/pi/(d/12)^2;
    CDll(i) = CDll(i)+dCDpro(i);

    if mach(i)<0.78
        Ke(i) = 0.00038;
    elseif mach(i)>=0.78 && mach(i)<=1.04
        Ke(i) = -0.4501*mach(i)^4+1.5954*mach(i)^3-2.1062*mach(i)^2+1.2288*mach(i)-0.26717;
    elseif mach(i)>1.04
        Ke(i) = 0.0002*mach(i)^2-0.0012*mach(i)+0.0018;
    else
%         error('mach number fail');
    end
    CDfe(i) = Ke(i)*4*((SB+Afe2)+4*S*(Cr+Ct))/pi/(d/12)^2;
    CDfe(i) = 0;
%     CDfe(i) = 0.5*CDfe(i);

    CDf(i) = CDfb(i)+1.04*CDff(i)+1.04*CDfpro(i)+CDfe(i)+1.04*CDll(i);
%     CDf(i) = CDf(i)/beta;
    
    if mach(i) <= 0.6 && t<tb
       CDf6 = CDfb(i)+1.04*CDff(i)+1.04*CDfpro(i)+CDfe(i)+1.04*CDll(i);
       CDb6 = Kb*((dr/(d/12))^nb)/sqrt(CDf6);
    end
    
    if mach(i)<=0.6
        CDb(i) = Kb*((dr/(d/12))^nb)/sqrt(CDf(i));
    elseif mach(i)>0.6
        if mach(i)<=1
            fb(i) = 1+215.8*(mach(i)-0.6)^6;
        elseif mach(i)<=2
            fb(i) = 2.0881*(mach(i)-1)^3-3.7938*(mach(i)-1)^2+1.4618*(mach(i)-1)+1.883917;
        elseif mach(i)>2
            fb(i) = 0.297*(mach(i)-1)^3-0.7937*(mach(i)-1)^2-0.1115*(mach(i)-1)+1.64006;
        end
        CDb(i) = CDb6*fb(i);
    end
    if t<=tb
%         CDb(i) = CDb(i)*0.8;
%         CDb(i) = 0;
    end

%     CDLEf(i) = n*S*ft/Aref*(((1+(gamma-1)/2*(mach(i)*cos(thetaLE))^2)^(gamma/(gamma-1))-1)...
%         /(gamma/2*(mach(i)*cos(thetaLE))^2));  %Too high
    
    CD(i) = CDf(i)+CDb(i)+CDt(i)+CDs(i);
        
%     CD0(i,1) = (1+60/(Rl/(d/12))^3+0.0025*LB/(d/12))*...
%                (2.7*Ln/(d/12)+4*LB/(d/12)+2*(1-dr/(d/12))*Lt/(d/12))*Cf(i,1);
%            
%     CD0(i,2) = 0.029*(dr/(d/12))^3/sqrt(CD0(i,1));
%     
%     CD0(i,3) = 2*Cf(i,2)*(1+2*Tf/Lf)*4*n*Afp/pi/dfa^2;
%     
%     CD0(i,4) = 2*Cf(i,2)*(1+2*Tf/Lf)*4*n*(Afp-Afe)/pi/dfa^2;
    
%     if x(i,2)==h && v(i,2)==0
%         phi(i) = 0;
%     else
%         phi(i) = asin((vwg(i,1)-v(i,1))/norm(vrw(i,:)));
%         phirot(i,:) = asin((vwg(i,1)-v(i,1)-vrotx(i,:))./vrwrot(i,:));
%     end
% 
%     aoa(i) = phi(i) - p(i);
%     aoai(i,:) = phirot(i,:) - p(i);
    
    if abs(aoa(i))<20*pi/180 && abs(aoa(i))>=4*pi/180
        aoa1 = 2*floor((abs(aoa(i))*180/pi)/2);
        aoa2 = aoa1+2;
        i1 = aoa1/2-1;
        i2 = i1+1;
        del = aoadd(i2)-(aoadd(i2)-aoadd(i1))*(aoa2-abs(aoa(i))*180/pi)/(aoa2-aoa1);
        ne = aoadn(i2)-(aoadn(i2)-aoadn(i1))*(aoa2-abs(aoa(i))*180/pi)/(aoa2-aoa1);
    elseif abs(aoa(i))<=4*pi/180
%         del = aoadd(1)+(aoadd(2)-aoadd(1))*(abs(aoa(i))*180/pi)/2;
%         ne = aoadn(1)+(aoadn(2)-aoadn(1))*(abs(aoa(i))*180/pi)/2;
        del = aoadd(1);
        ne = aoadn(1);
    else
        del = aoadd(9);
        ne = aoadn(9);
    end
    
    if abs(aoa(i))<=20*pi/180 %20*pi/180
        CDab = 2*del*(aoa(i))^2+3.6*ne*(1.36*Rl-0.55*Ln)*abs(aoa(i))^3/(pi*(d/12));
        CDaf = (aoa(i))^2*(1.2*Afp*4/(pi*dfa^2)+3.12*(kfb+kbf-1)*Afe*4/(pi*dfa^2));
    else
        CDab = 2*del*(20*pi/180)^2+3.6*ne*(1.36*Rl-0.55*Ln)*(20*pi/180)^3/(pi*(d/12));
        CDaf = (20*pi/180)^2*(1.2*Afp*4/(pi*dfa^2)+3.12*(kfb+kbf-1)*Afe*4/(pi*dfa^2));
    end
    if dist(i) < RL
        CDab = 0;
        CDaf = 0;
    end
%     CDab = 0;
%     CDaf = 0;
    CDabs(i) = CDab;
    CDafs(i) = CDaf;
    
    cl8(i) = sum(CNa(i,:))*aoa(i)*Ecl;
    %cl(i) = CNan*aoan(i)+CNat*aoat(i)+CNaf*aoaf(i)+CNab*aoab(i);
%     cl8(i) = sum(CNa(i,:).*aoai(i,:));
    %cd8(i) = sum(CD0(i,:))+CDab+CDaf;
    
%     CDfd(i) = (Afe/Aref)*cl8(i)^2/(pi*Afin);
%     CDfd(i) = CDfd(i)/beta;
%     if mach(i) < 1.0
%         e = 0.9;
%     else
%         e = 0.6;
%     end
        
%     CDa(i) = cl8(i)^2/(pi*e*Afin)/beta;
    CDa(i) = 5*sin(abs(aoa(i)));
    CDa(i) = 0;
    cd8(i) = (CD(i)+CDa(i))*Ecd;
    cd(i) = cd8(i);
%     cd8(i) = CD(i)+(CDabs(i)+CDafs(i))/beta;
%     cd8(i) = cd(i)+CDab+CDaf;
    %cd8(i) = CD(i)+CDfd(i);
%     cd8(i) = cd8(i)*cos(abs(aoa(i)));
    
%     cd8(i) = (cd8(i)*cos(abs(aoa(i)))-0.5*cl8(i)*sin(2*abs(aoa(i))))/(1-sin(abs(aoa(i)))^2);
%     cd8(i) = (cd8(i)-cl8(i)*sin(abs(aoa(i))))/cos(abs(aoa(i)));
    if cd8(i) < 0.3
        cd8(i) = 0.3;
    elseif cd8(i) > 1.6
        cd8(i) = 1.6;
    end
%     cl8(i) = cl8(i)/beta+0.1219;
%     cd8(i) = cd8(i)/beta-0.1228;
%     cl8(i) = cl8(i)/beta;
    %cd8(i) = cd8(i)/beta;
    if Rnc(i,1) == 0
        cd8(i) = 0;
    end
    
%     Xf = CPi(i,4);
% 
%     rf = (r0/12+DR)/2+S/3;
%     
%     for fi=1:1:4
%        PiB(fi,:) = [rf*sin(-pi/2+fi*pi/2),rf*cos(-pi/2+fi*pi/2),-Xf];
%        Pi(fi,:) = transpose(RM*transpose(PiB(fi,:)))+x(i,:);
%        Si(fi,:) = x(i,:)-Pi(fi,:);
% %        Si(fi,:) = transpose(-RM*transpose(PiB(fi,:)));
%        if norm(omega(i,:))~=0
%            Vp(fi,:) = norm(omega(i,:))*norm(Si(fi,:))*sin(acos(dot((Si(fi,:)/norm(Si(fi,:))),...
%                (omega(i,:)/norm(omega(i,:))))))*cross((Si(fi,:)/norm(Si(fi,:))),...
%                (omega(i,:)/norm(omega(i,:))))/norm(cross((Si(fi,:)/norm(Si(fi,:))),...
%                (omega(i,:)/norm(omega(i,:)))))+vrw(i,:);
%        else
%            Vp(fi,:) = vrw(i,:);
%        end
% %        li(fi,:) = [cos(-pi/2+fi*pi/2),sin(-pi/2+fi*pi/2),0];
%        if fi==1
%            li(fi,:) = transpose(pA(:,i));
%        elseif fi==2
%            li(fi,:) = transpose(yA(:,i));
%        elseif fi==3
%            li(fi,:) = transpose(-pA(:,i));
%        elseif fi==4
%            li(fi,:) = transpose(-yA(:,i));
%        end
%        Qc(1,i) = cos(cant/2);
%        Qc(2:4,i) = sin(cant/2)*cross(li(fi,:),rA(:,i))/norm(cross(li(fi,:),rA(:,i)));
%        Rc = [1-2*Qc(3,i)^2-2*Qc(4,i)^2, 2*Qc(2,i)*Qc(3,i)-2*Qc(1,i)*Qc(4,i),...
%            2*Qc(2,i)*Qc(4,i)+2*Qc(1,i)*Qc(3,i);...
%            2*Qc(2,i)*Qc(3,i)+2*Qc(1,i)*Qc(4,i),...
%            1-2*Qc(2,i)^2-2*Qc(4,i)^2, 2*Qc(3,i)*Qc(4,i)-2*Qc(1,i)*Qc(2,i);...
%            2*Qc(2,i)*Qc(4,i)-2*Qc(1,i)*Qc(3,i),...
%            2*Qc(3,i)*Qc(4,i)+2*Qc(1,i)*Qc(2,i), 1-2*Qc(2,i)^2-2*Qc(3,i)^2];
%        aoafi(i,fi) = pi/2-acos(dot((Vp(fi,:)/norm(Vp(fi,:))),Rc*transpose(li(fi,:))));
%        Vpn(i,fi) = norm(Vp(fi,:));
%     end
%     aoaf(i) = sum(aoafi(i,:))/n;
%     Vf(i) = sum(Vpn(i,:))/n;
% %     CR(i) = aoaf(i)/750;
%     LFr(i) = (CNaf/n)*aoaf(i);
%     
%     rf = (r0/12+DR)/2;
%     CNa0 = 2*pi/beta;
%     Clf(i) = n*(MAC+rf)*(CNaf/n)*cant/(d/12);
%     term = (Cr+Ct)/2*rf^2*S+(Cr+2*Ct)/3*rf*S^2+(Cr+3*Ct)/12*S^3;
%     if mach(i)<1
%         rollrate(i) = Aref*beta*norm(vrw(i,:))*MAC*(CNaf/n)*cant/(2*pi*term);
%     elseif mach(i)>=1
%         rollrate(i) = Aref*beta*norm(vrw(i,:))*MAC*(CNaf/n)*cant/(2*pi*term);
%     end
%     Cdf(i) = n*CNa0*rollrate(i)*term/Aref/(d/12)/norm(vrw(i,:));
%     CR(i) = Clf(i)-Cdf(i);
%     rf = (r0/12+DR)/2+S/3;
    
    
    %L(i,:) = cl8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*[cos(p(i)),sin(p(i))];
    
%     Lx(i) = L(i)*cos(p(i));
%     Lz(i) = L(i)*sin(p(i));
    
%     D(i,:) = cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*[sin(phi(i)),-cos(phi(i))];
% %     Dx(i) = D(i)*sin(phi(i));
% %     Dz(i) = -D(i)*cos(phi(1));
    if v(i,3) > 0
%         D(i,:) = cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*[sin(phi(i)),-cos(phi(i))];
        D(i,:) = -cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*rA(:,i);
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
%             cd8(i) = 1.6;
%             d = 24;
            cd8(i) = 2.2;
            d = 60;
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
%             cd8(i) = 1.6;
%             d = 24;

            Aref = pi*(d/12/2)^2;
%             D(i,:) = cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*[sin(phi(i)),-cos(phi(i))];
            D(i,:) = -cd8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*vhat(i,:);
            %             D(i) = cd(i)/2*A*rho(i)*(vrw(i))^2;
            %             Dx(i) = D(i)*sin(phi(i));
            %             Dz(i) = -D(i)*cos(phi(i));
        end
    end
    
    L(i,:) = cl8(i)/2*Aref*rho(i)*(norm(vrw(i,:)))^2*...
        cross(rA(:,i),cross(rA(:,i),vrw(i,:)/norm(vrw(i,:))))/...
        norm(cross(rA(:,i),cross(rA(:,i),vrw(i,:)/norm(vrw(i,:)))));
    Q(i) = 0.5*rho(i)*(norm(vrw(i,:)))^2;
    Qa(i) = Q(i)*aoa(i);
    QCNa(i) = Q(i)*sum(abs(CNa(i,:)));
    if v(i,3) < 0
        L(i,:) = [0,0,0];
    end    

    if t<=tb
        W(i,:) = [0,0,-w(i)];
    else
        W(i,:) = [0,0,-dw];
    end
    Eomega1 = [0,Eomega,0];
    C(i,:) = 2*m*cross(Eomega1,v(i,:))-m*cross(Eomega1,cross(Eomega1,x(i,:)));
    %C(i,:) = [0,0];
    
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
    
    if Iapogee2 ~= 0 && ejection==0 && cd8(i) ~= 2.2
        Shock = [1700/3,1700];
        Reac(i,:) = Shock(1)*[-sign(v(i,1))*cos(p(i)),0,-sin(p(i))]; %Make perpendicular to axis
        Reac(i,:) = Reac(i,:)+Shock(2)*[-sin(p(i)),0,cos(p(i))]; %Axial
        Mreac = [0,-sign(p(i))*Shock(1)*(cg(i)-47/12),0];
        ejection = 1;
        % [-sin(p(i)),cos(p(i))]
    else
        Mreac = [0,0,0];
    end
    
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
%     aA(i) = AF(i)/m;
%     aN(i) = NF(i)/m;

    cgp(i) = (cgf(i)*wf(i)+cgo(i)*wo(i))/(wf(i)+wo(i));
    Mt(i,:) = -(wfr(i)/g)*((Rl-cg(i))^2-(cgp(i)-cg(i))^2)*...
        transpose(RM*diag([1 1 0],0)/RM*transpose(omega(i,:)));
    Mt(i)=0;
    
    if CP8(i)<cg(i) && dist(i)>RL-Rl+cg(i) && v(i,3)>0 && abs(aoa(i))<15*pi/180
%         disp(i);
        %error('Unstable');
    end

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
    MT(i,:) = -(Rl-cg(i))*cross(rA(:,i),Tr(i,:));
%     MT(i,:) = [0,0,0];
%     M(i) = MN(i)+Mt(i)+Md(i)+Mll(i)+Mpro(i)+Mreac;
%     Mr(i,:) = CR(i)/2*rho(i)*Aref*(rollrate(i)*rf)^2*(d/12)*transpose(rA(:,i));
%     Mr(i,:) = -CR(i)*LFr(i)*rf*transpose(rA(:,i));
    Mr(i,:) = [0,0,0];
    M(i,:) = MN(i,:)+Mt(i,:)+Mll(i,:)+Mpro(i,:)+Mreac+MT(i,:)+Mr(i,:);
%     alpha(i,:)=(I0)\transpose(M(i,:)-transpose(Idot*transpose(omega(i,:))));
    alpha(i,:)=(I0)\(transpose(M(i,:)-...
        cross(omega(i,:),I0*transpose(omega(i,:)))));
    if dist(i)<RL || v(i,3)<0
        alpha(i,:) = [0,0,0];
    end
    omega(i+1,:) = omega(i,:) + alpha(i,:)*dt;
    y(i+1) = y(i) + omega(i,1)*dt;
    p(i+1) = p(i) + omega(i,2)*dt;
    r(i+1) = r(i) + omega(i,3)*dt;
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
    gs(i) = g;
    ms(i) = m;
%     if dist(i) > RL || v(i,2) < 0 && v(i,2)>0
%         NF2(i) = sign(alpha(i))*norm(cross(F(i,:),[-sin(p(i)),cos(p(i))]));
%     end

    iters(i) = i;
    i = i + 1;
end
xs(si,1:3) = x(itf,:);
[apogee,Iapogee] = max(zs);
apogees(si) = zs(Iapogee);
OTRSs(si) = OTRs;
deployment_vs(si) = s(Iapogee);
drift(si) = norm(x(itf,:)-x(1,:));
maxS(si) = max(s);
end
xs(:,1:2) = -xs(:,1:2)/5280;
x(:,1:2) = -x(:,1:2)/5280;
hold on
plot(xs(2:end,1),xs(2:end,2),'x','MarkerSize',10,'LineWidth',2.0);
hold on
plot(x(1,1),x(1,2),'o','MarkerSize',20,'LineWidth',2.0);
hold on
plot(xs(1,1),xs(1,2),'+','MarkerSize',20,'LineWidth',2.0);

%% Obtain Statistical Landing Plots

hold on
for scp = 0.95:-0.05:0.8
    sc = -2*log(1-scp);
    [evec,eval] = eig(sc*cov([xs(2:end,1)-xs(1,1) xs(2:end,2)-xs(1,2)]));
    if eval(1,1)<=eval(2,2)
        seval = eval(1,1);
        sevec = evec(1,:);
        leval = eval(2,2);
        levec = evec(2,:);
    else
        seval = eval(2,2);
        sevec = evec(2,:);
        leval = eval(1,1);
        levec = evec(1,:);
    end
    angle = atan2(levec(2),levec(1));
    if angle<0
        angle = angle+2*pi;
    end
    angles = linspace(0,2*pi);
    ellipse = transpose([sqrt(sc*leval)*cos(angles);sqrt(sc*seval)*sin(angles)])*...
        [cos(angle),sin(angle);-sin(angle),cos(angle)];
    plot(ellipse(:,1)+xs(1,1),ellipse(:,2)+xs(1,2),'-','LineWidth',2.0);
end

xlabel('Easting (mi)','FontSize',14);
ylabel('Southing (mi)','FontSize',14);
title('Landing Positions','FontSize',14);
legend({'Stochastic Landings','Launch Pad','Mean Landing','95% Confidence',...
    '90% Confidence','85% Confidence','80% Confidence'},'FontSize',14);


