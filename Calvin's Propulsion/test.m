clear all; clc;

addpath('Classes')
addpath(genpath('Supporting Functions'))
load('95EthLOxCompact')


He = Gas;
He.R = 12421; %[ft lbf/slug °R]
He.gam = 1.66;
He.G = 0.138;
He.c_v = He.R/(He.gam - 1);
He.c_p = He.R*He.gam/(He.gam - 1);

LOx = Liquid;
LOx.rho = 71.17/32.2;

Eth = Liquid;
Eth.rho = 50.11/32.2;

CV_o = CheckValve;
CV_o.name = "CV_o";
CV_o.Cv = 1.7;
CV_o.p_crack = 5;
CV_o.r_pop = 0.5;

CV_f = CheckValve;
CV_f.name = "CV_f";
CV_f.Cv = 1.7;
CV_f.p_crack = 5;
CV_f.r_pop = 0.5;

Reg = Regulator;
Reg.name = "Reg";
Reg.Cv = 0.2;
Reg.p_set = 500;
Reg.p_crack = 25;
Reg.r_droop = (450-325)/(20-140);
Reg.p_droop_cali = 2400;

PT = PressurantTank;
PT.name = "PT";
PT.p = 4500;
PT.T = 519;
PT.V = 100/12^3;
PT.rho = IdealGasDensity(PT.p,PT.T,He.R);
PT.m = PT.rho*PT.V;

IntT = Tubing;
IntT.name = "IntT";
IntT.p = 500;
IntT.T = 519;
IntT.V = 1/12^3;
IntT.rho = IdealGasDensity(PT.p,PT.T,He.R);
IntT.m = PT.rho*PT.V;

OxT = PropellantTank;
OxT.name = "OxT";
OxT.p = 500;
OxT.T = 519;
OxT.V = 1/12^3;
OxT.rho = IdealGasDensity(OxT.p,OxT.T,He.R);
OxT.m = OxT.rho*OxT.V;
OxT.m_prop = 0.833;
OxT.rho_prop = LOx.rho;

FuT = PropellantTank;
FuT.name = "FuT";
FuT.p = 500;
FuT.T = 519;
FuT.V = 1/12^3;
FuT.rho = IdealGasDensity(FuT.p,FuT.T,He.R);
FuT.m = FuT.rho*FuT.V;
FuT.m_prop = 0.965;
FuT.rho_prop = Eth.rho;

Inj = InjectorManifold;
Inj.name = "OxMan";

map = {};
map = MapInSeries(map,PT);
map = MapInSeries(map,Reg,PT);
map = MapInSeries(map,IntT,Reg);
%map = MapInSeries(map,CV_o,IntT);
map = MapInSeries(map,CV_o,IntT);
map = MapInParallel(map,CV_f,CV_o);
map = MapInSeries(map,OxT,CV_o);
map = MapInSeries(map,FuT,CV_f);
map = MapInSeries(map,Inj,OxT);
map = MapInSeries(map,Inj,FuT);

% map = MapInParallel(map,CV_o,Reg);
% map = MapInSeries(map,FuT,CV_o);

%map = DefineInletsOutletsEquations(map);

%DefineGoverningEquations(map,He);

volumes = VectorizeMapVolumes(map);

tf = 20;
dt = 0.001;
nsim = ceil(tf/dt);
p_PT = zeros(1,nsim);
T_PT = zeros(1,nsim);
m_PT = zeros(1,nsim);
p_IntT = zeros(1,nsim);
T_IntT = zeros(1,nsim);
m_IntT = zeros(1,nsim);
p_OxT = zeros(1,nsim);
Tul_OxT = zeros(1,nsim);
mul_OxT = zeros(1,nsim);
p_FuT = zeros(1,nsim);
Tul_FuT = zeros(1,nsim);
mul_FuT = zeros(1,nsim);

i = 0;
t = 0;
while (OxT.m_prop > 0 && FuT.m_prop > 0)
    fprintf('i = %d\n',i)
    i = i + 1;
    p_PT(i) = PT.p;
    T_PT(i) = PT.T;
    V_PT(i) = PT.V;
    m_PT(i) = PT.m;
    
    p_IntT(i) = IntT.p;
    T_IntT(i) = IntT.T;
    V_IntT(i) = IntT.V;
    m_IntT(i) = IntT.m;
    
    p_OxT(i) = OxT.p;
    Tul_OxT(i) = OxT.T;
    Vul_OxT(i) = OxT.V;
    mul_OxT(i) = OxT.m;
    
    p_FuT(i) = FuT.p;
    Tul_FuT(i) = FuT.T;
    Vul_FuT(i) = FuT.V;
    mul_FuT(i) = FuT.m;
    
    fprintf('%.2f\n',PT.p)
    fprintf('%.2f\n',IntT.p)
    fprintf('%.2f\n',OxT.p)
    fprintf('%.2f\n',FuT.p)
    
    Xvals = VectorizeXvals(volumes);
    
    DefineGoverningEquations(map,He);
    
    A = BlockAMatrix(volumes);

    mdots = MassFlowRates(map,Xvals,A,dt,He,CombData);
    bvals = mdots;
    
    dX = A\bvals;
    
    Xvals = Xvals + dX*dt;
    
    DealXvals(volumes,Xvals);
    
    if isreal(Xvals)
        fprintf('real\n')
    else
        error('values complex')
    end
    
    t = t + dt;
end
    
lw = 2;
figure
hold on
plot(p_PT(1:i),'LineWidth',lw)
plot(p_IntT(1:i),'LineWidth',lw)
plot(p_OxT(1:i),'LineWidth',lw)
plot(p_FuT(1:i),'LineWidth',lw)
hold off
legend('PT','IntT','OxT','FuT')
title('pressure')

figure
hold on
plot(T_PT(1:i),'LineWidth',lw)
plot(T_IntT(1:i),'LineWidth',lw)
plot(Tul_OxT(1:i),'LineWidth',lw)
plot(Tul_FuT(1:i),'LineWidth',lw)
hold off
legend('PT','IntT','OxT','FuT')
title('temperature')

figure
hold on
plot(V_PT(1:i),'LineWidth',lw)
plot(V_IntT(1:i),'LineWidth',lw)
plot(Vul_OxT(1:i),'LineWidth',lw)
plot(Vul_FuT(1:i),'LineWidth',lw)
hold off
legend('PT','IntT','OxT','FuT')
title('volume')

figure
hold on
plot(m_PT(1:i),'LineWidth',lw)
plot(m_IntT(1:i),'LineWidth',lw)
plot(mul_OxT(1:i),'LineWidth',lw)
plot(mul_FuT(1:i),'LineWidth',lw)
plot(m_PT+m_IntT+mul_OxT+mul_FuT,'LineWidth',lw)
hold off
legend('PT','IntT','OxT','FuT','Total')
title('mass')