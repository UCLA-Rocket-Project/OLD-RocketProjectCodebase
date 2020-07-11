function map = AresSystemMapping(fluid)

He = fluid.He;
LOx = fluid.LOx;
Eth = fluid.Eth;

%% Valves and Regulator
Reg = Regulator;
Reg.name = "Reg";
Reg.Cv = 0.2;
Reg.p_set = 450;
Reg.p_crack = 25;
Reg.r_droop = (450-325)/(20-140);
Reg.p_droop_cali = 2400;
Reg.SPE = 0.025;

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

%% Volumes
PT = PressurantTank;
PT.name = "PT";
PT.p = 4500;
PT.T = 519;
PT.V = 120/12^3;
PT.rho = IdealGasDensity(PT.p,PT.T,He.R);
PT.m = PT.rho*PT.V;

IntT = Tubing;
IntT.name = "IntT";
IntT.p = Reg.p_set - Reg.SPE*(PT.p - Reg.p_set);
IntT.T = 519;
IntT.V = 5/12^3;
IntT.rho = IdealGasDensity(PT.p,PT.T,He.R);
IntT.m = PT.rho*PT.V;

OxT = PropellantTank;
OxT.name = "OxT";
OxT.p = Reg.p_set - Reg.SPE*(PT.p - Reg.p_set);
OxT.T = 519;
OxT.V = 10/12^3;
OxT.rho = IdealGasDensity(OxT.p,OxT.T,He.R);
OxT.m = OxT.rho*OxT.V;
OxT.m_prop = 0.833;
OxT.rho_prop = LOx.rho;

FuT = PropellantTank;
FuT.name = "FuT";
FuT.p = Reg.p_set - Reg.SPE*(PT.p - Reg.p_set);
FuT.T = 519;
FuT.V = 10/12^3;
FuT.rho = IdealGasDensity(FuT.p,FuT.T,He.R);
FuT.m = FuT.rho*FuT.V;
FuT.m_prop = 0.559;
FuT.rho_prop = Eth.rho;

%% Injector Manifolds
FuMan = InjectorManifold;
FuMan.name = "FuMan";
FuMan.type = "fuel";
FuMan.CdA = 0.02895139521*1.0;

OxMan = InjectorManifold;
OxMan.name = "OxMan";
OxMan.type = "oxidizer";
OxMan.CdA = 0.03041978109*1.2;

%% Combustion Chamber
CC = CombustionChamber;
CC.name = "CC";
CC.OF = OxMan.CdA*sqrt(OxT.rho_prop)/(FuMan.CdA*sqrt(FuT.rho_prop));
CC.A_t = pi*(1.5/2)^2;
CC.cT = 1.4;

%% Mapping
map = {};
map = MapInSeries(map,PT);
map = MapInSeries(map,Reg,PT);
map = MapInSeries(map,IntT,Reg);
map = MapInSeries(map,CV_o,IntT);
map = MapInParallel(map,CV_f,CV_o);
map = MapInSeries(map,OxT,CV_o);
map = MapInSeries(map,FuT,CV_f);
map = MapInSeries(map,OxMan,OxT);
map = MapInSeries(map,FuMan,FuT);
map = MapInSeries(map,CC,OxMan);
map = MapInSeries(map,CC,FuMan);

