%% Input Parameters

recovery = 0; % 0 for no recovery / 1 for recovery
turb = 1; % 0 for none / 1 for light / 2 for medium / 3 for severe

h = 3500;
x0 = 282;
y0 = 30.5;
Temp0 = 70;
kTemp = 10000; % Seems to be the most reasonable "tending to atmospheric model" rate
P0 = 14.24*144;
hr0 = 0;
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
vwgx = -5*5280/3600;
vwgy = -5*5280/3600;
k = 0.4;
dref = 0.000;
z1 = 0.1/100;
atm_stability = 1; % 1=neutral, 2=stable, 3=unstable

dw = 91.03;
T0 = readmatrix('thrustArray.xlsx');
% T0 = readmatrix('Odyssey3Thrust.xlsx');
% T0 = readmatrix('PropSimThrust.xlsx');
Imp = trapz(T0(:,1),T0(:,2));
Impfactor = 9200/Imp-0.037;
T0 = T0*Impfactor;
Imp = trapz(T0(:,1),T0(:,2));
% Isp = 200;
wp = 50;
OF = 1.5;
wfi = wp/(OF+1);
woi = wfi*OF;
ww = dw+wfi+woi;
Isp = Imp/(ww-dw);
cgfi = 124.43/12;
cgoi = 96.18/12;
lfi = 22.25/12;
loi = 25.07/12;
xfuel0 = cgfi-lfi/2;
xox0 = cgoi-loi/2;
dtank = 7/12;
rtank = dtank/2;
% dnoz = 4;
dnoz = 2.81;
cgd = 101.28/12;
RL = 55;
RLa = 3*pi/180;
p0 = RLa;
RLv = [-RL*sin(RLa),RL*cos(RLa)];
df = (wfi/g)/(lfi*pi*rtank^2);
dox = (woi/g)/(loi*pi*rtank^2);
ejection = 0;

ti = 0;
tb = T0(end,1);
% tb = round(find(T0(:,2),1,'last')/1000,2); %PropSim
tf = 350;
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
Ln = 35/12;
xn = 0.466*Ln;
LB = 125.23/12;
xp = Ln+LB;
Lt = 15.75/12;
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
% Apro = 1/144;
Apro = 1.47/144;
Lp = 57/12;
aPro = 70.32;
rpro = sqrt(Apro/pi);
Spro = 2*2*Apro+2*pi*rpro*Lp;
% Spro = 2*2*Apro+2*2*pi*rpro/2*Lp;
% Spro = 113.03/144;

All = 1/144; %in
% All = 0/144; %in
rll = sqrt(All/pi);
Lll = 6/12;
Sll = 2*All+2*2*pi*rll*Lll;
aLL = Rl/2; %Midway on rocket
muk = 0.12; %Steel on steel lubricated
mus = 0.16; %Steel on steel lubricated
% muk = 0; %Steel on steel lubricated
% mus = 0; %Steel on steel lubricated
if p0<0
%     mus = -mus;
end

Lc = [Rl,Cr,Lp,Lll];

%Find cg formulas later
cgn = 24.36/12;
cgb = Ln+LB/2;
cgt = 166/12;
cgfin = 172/12;
cgi = [cgn, cgb, cgt, cgfin];

vwg = zeros(t_steps,3); % Wind model in loop
ustar = zeros(t_steps,1);
vstar = zeros(t_steps,1);
% ustar(:) = abs(vwgx)/3.28*k/log((RLv(2)/3.28-dref/3.28)/(z1/3.28));
cv = 0.25; % 0.1 to 0.4
fv = 2*Eomega*sin(lat);
% ustar(:) = vwgx*k/log(RLv(2)/z1);
if atm_stability == 1
    Lv = 500*3.281;
    ustar(:) = abs(vwgx)*k/log(RLv(2)/z1);
    vstar(:) = abs(vwgy)*k/log(RLv(2)/z1);
    hv = cv*ustar(1)/abs(fv);
    hv2 = cv*vstar(1)/abs(fv);
    vA = 4.9;
    vB = 1.9;
    LMBL = hv/2*(((log(ustar(1)/abs(fv)/z1)-vB)^2+vA^2)^0.5-log(hv/z1))^(-1);
    LMBL2 = hv2/2*(((log(vstar(1)/abs(fv)/z1)-vB)^2+vA^2)^0.5-log(hv2/z1))^(-1);
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
if vwgx<0
    ustar(:) = -ustar(:);
end
if vwgy<0
    vstar(:) = -vstar(:);
end

% tgust = 5.82;
tgust = 15.93;
dtgust = 3;
dtrise = 0*dtgust;
gust = -30;
% gust = 0;

% tma = [0.16/180*pi,0];
% tmap = 0.0;
tmap = 0.7;
tmam = 0.32/180*pi;
tma = tmam*[cos(tmap),sin(tmap)];
% tma = 0*[sin(tmap),cos(tmap)];

% cant = 1/180*pi;
cant = 0;
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
x(1,1) = x0;
x(1,2) = y0;
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
Id = 1239.3;
Idz = 3.9;
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

%% Save Input File

save('InputParams_6DOF');
