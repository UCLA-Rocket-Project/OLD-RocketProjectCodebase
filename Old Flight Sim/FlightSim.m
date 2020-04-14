%Flight Sim (Problem 3 + 4 component form)
%last edited 10/21/2019
%Problem 3
%Units:sec,in,slugs,lbf
clear variables; clc; close all;
%% General Information
    g = 32.1740;
    k = 1.5;
        %coefficient in body lift force equation, between 1.1 and 1.5
    t0 = 0; tf = 60; dt = 0.01;
    [T_data, tb] = getThrust('2-22-20.xlsx', dt); % need thrust data for tb unless we assume ideal tb
    steps = ((tf-t0)/dt)+1; tbsteps = ((tb-t0)/dt)+1;
    t = linspace(t0,tf,steps); % time (s)
        
    %% Parachutes
    %Drogue
    drogueDiameter = 24/12; %unit: FT
    drogueCD = 1.6;
    drogueArea = pi*(drogueDiameter/2)^2;
    mainDiameter = 60/12; %unit: FT
    mainCD = 2.2;
    mainArea = pi*(mainDiameter/2)^2;
    mainAlt = 1500; %deployment altitude
    
    %% Ogive Nose
    noseDiameter = 7;
    FinenessRatio = 5;
    noseLength = noseDiameter*FinenessRatio;
    noseWeight = 2.74; %
    noseMass = noseWeight/g;
    noseThickness = 0.3;
    noseCG = 24.36;
    noseAP = 0.5*noseLength*noseDiameter;
        %modeled as a conical shape
    
    %% Payload (Recovery, avionics, electronics), assume d = bodyDiameter
    payloadWeight = 15.6574; %includes recovery, payload, electronics, associated housing and bulkheads
    payloadMass = payloadWeight/g;
    payloadLength = 20.5;
    payloadX = noseLength;
    payloadCG = 39.5595; % caluclated from "2019-2020 Comprehensive CG Calculator,"
                       % ~ payloadX + payloadLength/2
   
    
    %% Body
    bodyDiameter = 7;
    bodyLength = 127.7;
    bodyWeight = 21.7825; %includes misc. prop, longerons, aeroshell, bolts etc, press bulkhead, rail guides, couplers, thrust bulkhead, bulkhead brackets
    bodyMass = payloadWeight/g;
    bodyThickness = 0.05;
    bodyKfactor = 1.5;
    bodyCG = (bodyLength/2)+noseLength;
    bodyAP = bodyLength*bodyDiameter;
    
    %% Pressurant Tank
        %assume diameter = bodyDiameter, no assignment
    presstWeight = 2.25;
    presstMass = presstWeight/g;
    presstLength = 9.5;
    presstThickness = 1;
    presstX = 55.5;
        %"X" is distance from top of component to nose tip
    presstCG = (presstLength/2)+presstX;
    
    %% Plumbing (modeled as cylinders, d = bodyDiameter)
    %Press-Ox Plumbing
    poplumbWeight = 6.59;
    poplumbMass = poplumbWeight/g;
    poplumbLength = 10;
    poplumbX = 65;
    poplumbCG = (poplumbLength/2)+poplumbX;
    %Ox-Fuel Plumbing
    ofplumbWeight = 1.49;
    ofplumbMass = ofplumbWeight/g;
    ofplumbLength = 5.5;
    ofplumbX = 106.5;
    ofplumbCG = (ofplumbLength/2)+ofplumbX;
    %Fuel-Engine Plumbing
    feplumbWeight = 4.35;
    feplumbMass = feplumbWeight/g;
    feplumbLength = 14;
    feplumbX = 134;
    feplumbCG = (feplumbLength/2)+feplumbX;
    
    %% Oxidizer Tank
        %assume diameter = bodyDiameter
    oxtWeight = 6.22;
    oxtMass = oxtWeight/g;
    oxtLength = 31.5;
    oxtThickness = 0.09;
    oxtcapThickness = 4.06;
    oxtX = 75;
    oxtCG = (oxtLength/2)+oxtX;
    
    %% Fuel Tank
        %assume diameter = bodyDiameter
    fueltWeight = 5.5;
    fueltMass = fueltWeight/g;
    fueltLength = 27.5;
    fueltThickness = 0.09;
    fueltcapThickness = 4.06;
    fueltX = 112;
    fueltCG = (fueltLength/2)+fueltX;
    
    %% Propellant
    oxWeighti = 28;
    oxWeight = Emptying(oxWeighti,0,tbsteps,steps);
    oxMass = oxWeight/g;
    oxLengthi = oxtLength-(2*oxtcapThickness);
    oxLength = Emptying(oxLengthi,0,tbsteps,steps);
    oxXi = oxtX+oxtcapThickness;
    oxX = Emptying(oxXi,oxXi+oxLengthi,tbsteps,steps);
    oxCG = oxX+(0.5*oxLength);
    fuelWeighti = 18.5;
    fuelWeight = Emptying(fuelWeighti,0,tbsteps,steps);
    fuelMass = fuelWeight/g;
    fuelLengthi = fueltLength-(2*fueltcapThickness);
    fuelLength = Emptying(fuelLengthi,0,tbsteps,steps);
    fuelXi = fueltX+fueltcapThickness;
    fuelX = Emptying(fuelXi,fuelXi+fuelLengthi,tbsteps,steps);
    fuelCG = fuelX+(0.5*fuelLength);    
    
    %% Engine
    engWeight = 14.79; %includes engine centering rings, engine
    engMass = engWeight/g;
    engLength = 19;
    engDiameter = 4; %%%%%%%
    engX = 153.5;
    engCG = (engLength/2)+engX;
    
    %% Fins (Clipped Delta)
    finRC = 14.5;
    finTC = 3.25;
    finSS = 8.5;
    finLESA = atan((finRC-finTC)/finSS);
        %Leading edge sweep angle (clipped delta formula);
    finMC = sqrt(((finRC-finTC)/2)^2 + finSS^2);
        %Mid-chord length
    finNumber = 4;
    finWeight = 0.9;
    finMass = finWeight/g;
    finsWeight = finWeight*finNumber;
    finsMass = finMass*finNumber;
    finsThickness = 0.7;
    finsX = 156.9;
    finsCG = finsX+(((2*finRC)/3) - ((finTC^2)/(3*(finRC+finTC))));
        %Clipped Delta derived formula from integration       
    
    %% Boattail
    tailLength = 9;
    tailDiameter = 5.4;
        r1 = bodyDiameter/2;
        r2 = tailDiameter/2;
    tailThickness = 0.2;
    tailWeight = 1.2;
    tailMass = tailWeight/g;
    tailX = noseLength+bodyLength;
    tailCG = tailX+((((r1/3)+((2*r2)/3)-(0.5*tailThickness))/(r1+r2-tailThickness))*tailLength);
        %Boattail derived formul from integration
    tailAP = 0.5*tailLength*(bodyDiameter+tailDiameter);
        
    %% General Rocket Mass Calcs
    DryWeight = noseWeight+bodyWeight+finsWeight+tailWeight...
                +payloadWeight+presstWeight+oxtWeight+fueltWeight+engWeight...
                +poplumbWeight+ofplumbWeight+feplumbWeight;
    DryMass = DryWeight/g;
    WetWeight = DryWeight+oxWeighti+fuelWeighti;
    WetMass = WetWeight/g;
    Wy = -(DryWeight+oxWeight+fuelWeight);
    m = -Wy/g;
    Length = noseLength+bodyLength+tailLength;
    CG=((noseMass*noseCG)+(bodyMass*bodyCG)+(tailMass*tailCG)+(finsMass*finsCG)...
        +(payloadMass*payloadCG)+(presstMass*presstCG)+(oxtMass*oxtCG)+(fueltMass*fueltCG)...
        +(engMass*engCG)+(oxMass.*oxCG)+(fuelMass.*fuelCG)...
        +(poplumbMass*poplumbCG)+(ofplumbMass*ofplumbCG)+(feplumbMass*feplumbCG))./m;
    AP = noseAP + bodyAP + tailAP;

%% Moment of Inertia Calculations
    nosePtMOI = PointMOI(noseMass,noseCG,CG);
    noseMOI = nosePtMOI + (noseMass/4)*((noseDiameter/2)^2 + 2*noseLength^2);
    
    payloadPtMOI = PointMOI(payloadMass,payloadCG,CG);
    payloadMOI = payloadPtMOI + CylinderMOI(payloadMass,bodyDiameter/2,0,payloadLength);
    
    bodyPtMOI = PointMOI(bodyMass,bodyCG,CG);
    bodyMOI = bodyPtMOI + CylinderMOI(bodyMass,bodyDiameter/2,bodyDiameter/2-bodyThickness,bodyLength);
    
    presstPtMOI = PointMOI(presstMass,presstCG,CG);
    presstMOI = presstPtMOI + CylinderMOI(presstMass,bodyDiameter/2,bodyDiameter/2-presstThickness,presstLength);
    
    poplumbPtMOI = PointMOI(poplumbMass,poplumbCG,CG);
    poplumbMOI = poplumbPtMOI + CylinderMOI(poplumbMass,bodyDiameter/2,0,poplumbLength);
    
    ofplumbPtMOI = PointMOI(ofplumbMass,ofplumbCG,CG);
    ofplumbMOI = ofplumbPtMOI + CylinderMOI(ofplumbMass,bodyDiameter/2,0,ofplumbLength);
    
    feplumbPtMOI = PointMOI(feplumbMass,feplumbCG,CG);
    feplumbMOI = feplumbPtMOI + CylinderMOI(feplumbMass,bodyDiameter/2,0,feplumbLength);
    
    oxtPtMOI = PointMOI(oxtMass,oxtCG,CG);
    oxtMOI = oxtPtMOI + CylinderMOI(oxtMass,bodyDiameter/2,bodyDiameter/2-oxtThickness,oxtLength);
    
    fueltPtMOI = PointMOI(fueltMass,fueltCG,CG);
    fueltMOI = fueltPtMOI + CylinderMOI(fueltMass,bodyDiameter/2,bodyDiameter/2-fueltThickness,fueltLength);
    
    oxPtMOI = PointMOI(oxMass,oxCG,CG);
    oxMOI = oxPtMOI + CylinderMOI(oxMass,bodyDiameter/2-oxtThickness,0,oxLength);
    
    fuelPtMOI = PointMOI(fuelMass,fuelCG,CG);
    fuelMOI = fuelPtMOI + CylinderMOI(fuelMass,bodyDiameter/2-fueltThickness,0,fuelLength);
    
    engPtMOI = PointMOI(engMass,engCG,CG);
    engMOI = engPtMOI + CylinderMOI(engMass,engDiameter/2,0,engLength);
    
    finsPtMOI = PointMOI(finsMass,finsCG,CG);
    %finsMOI = ClippedDeltaMOI(finRC,finTC,finSS,finWeight,finMass,finsX,bodyDiameter);
    %Suppressed until i figure out wtf is going on
    finsMOI = finsPtMOI;
    
    tailPtMOI = PointMOI(tailMass,tailCG,CG);
    %tailMOI = tailPtMOI + (tailMass/4)*((bodyDiameter/2)^2 + 2*tailLength^2);
    %Suppressed until i figure out wtf is going on
    tailMOI = tailPtMOI;
    
    PtMOI = nosePtMOI + bodyPtMOI + tailPtMOI + finsPtMOI... %% Not strictly necessary
          + payloadPtMOI + presstPtMOI + oxtPtMOI + fueltPtMOI... 
          + engPtMOI + oxPtMOI + fuelPtMOI...
          + poplumbPtMOI + ofplumbPtMOI + feplumbPtMOI;
    MOI = noseMOI + bodyMOI + tailMOI + finsMOI...
        + payloadMOI + presstMOI + oxtMOI + fueltMOI...
        + engMOI + oxMOI + fuelMOI...
        + poplumbMOI + ofplumbMOI + feplumbMOI;

%Normal force curve slope (summing up from aeroshell components)
    noseCNa = 2;
    tailCNa = 2*((tailDiameter/noseDiameter)^2 - 1);
    finsCNa = (1 + (bodyDiameter/2)/(finSS+bodyDiameter/2))...
             *((4*finNumber*(finSS/noseDiameter)^2)/...
             (1+sqrt(1+(2*finMC/(finRC+finTC))^2)));

%CP calculations
    noseCP = 0.466*noseLength;
    tailCP = tailX + tailLength*(1 + (1 - bodyDiameter/tailDiameter)/(1 - (bodyDiameter/tailDiameter)^2))/3;
    finCP = finsX + (finRC-finTC)*(finRC+2*finTC)/(3*(finRC+finTC))...
                  +(finRC + finTC - finRC*finTC/(finRC+finTC))/6;
    bodyCP = (1/AP)*(bodyAP*(noseLength + bodyLength/2) + 2*noseAP*noseLength/3 ...
                    + tailAP*(noseLength + bodyLength + (tailDiameter + 2*bodyDiameter)*tailLength/(3*(bodyDiameter+tailDiameter))));

CNai0 = [noseCNa; 0; tailCNa; finsCNa];
CPi = [noseCP; bodyCP; tailCP; finCP];
CGi = [noseCG; bodyCG; tailCG; finsCG];


% distMin = zeros(1,101);
% distMax = zeros(1,101);
% groundWind = zeros(1,101);
% ORaoa = zeros(1,101);
%for windex = 1:101
%% Problem 4 (component)
% General variables
RAS = readmatrix('CD Test.csv','Range','A2:E501');
    %Gets drag matrix
    Ma_dd = getDragDivergenceMach(RAS);
rail = 57 - (Length/12) + (CG(1)/12); %effective launch rail length(ft)
h0 = 2000; % (ft) FAR altitude

% Weather/Atmospheric Initializations (MonteCarlo start!)
%meanWind = %(windex-1)*(20*1.46667)/100;%20*1.46667;%(windex-1)*(20*1.46667/65.6533)/100; %wind velocity(mph to ft/s)
%groundWind(windex) = meanWind;%windProfileLog(wind,rail);
N_n = 50;
no_turbProb = 1; % Currently not using turbulence model
[mu_arr, sigma_arr, no_time] = getAtmosphereData('FlightSim Scratch.xlsm');
trials = no_time*no_turbProb*N_n;
[timeLaunch_arr, turbProb_arr, Temp0_arr, meanWind_arr, P0_arr] = setMonteCarlo(N_n, mu_arr, sigma_arr, no_time, no_turbProb); % set up initial condition arrays for N trials

% Initialize design variable arrays
    % Recovery Case-specific
    dist_arr = zeros(trials,1); %
    distT_arr = zeros(trials,1);
    
    % Can run for any recovery case
    tai_arr = zeros(trials,1);
    rx_arr = zeros(trials,steps);
    ry_arr = zeros(trials,steps); %
    vMax_arr = zeros(trials,1);
    MaMax_arr = zeros(trials,1);
    aMax_arr = zeros(trials,1);
    vDeploy_arr = zeros(trials,1); %
    vOR_arr = zeros(trials,1);
    aoaOR_arr = zeros(trials,1); %
    aoa_arr = zeros(trials,steps); %
    phi_arr = zeros(trials,steps);
for trialNo = 1:trials % Largest loop
    
% Set trial-specific variables
[C1, C2, Temp1] = setAtmosphere(h0,Temp0_arr(trialNo),P0_arr(trialNo));
meanWind = meanWind_arr(trialNo);
turbProb = turbProb_arr(trialNo);

% seed_ug = randi(2^32 - 1);
% seed_vg = randi(2^32 - 1);
% seed_wg = randi(2^32 - 1);
% seed_pg = randi(2^32 - 1);
% seed = [seed_ug, seed_vg, seed_wg, seed_pg];
%
% % Initialize Simulink Model
% set_param('windTurbulenceS/turbulenceGenerator', 'Wingspan', num2str((bodyDiameter + 2*finSS)/12));
% set_param('windTurbulenceS/turbulenceGenerator', 'W20', 'meanWind');                                                 
% set_param('windTurbulenceS/turbulenceGenerator', 'TurbProb', turbProb);
% set_param('windTurbulenceS/turbulenceGenerator','Seed','seed');

%Initial values
        A_e = 0.25*pi*(2.86/12)^2; %exit area of nozzle
    %Reference area(ft^2)
        A=pi*((bodyDiameter/12)/2)^2;

%Initialize arrays
aa = zeros(1,steps); %angular acceleration(rad/s^2)
    w = zeros(1,steps); %angular velocity(rad/s)
    p = zeros(1,steps); %angular position(rad) or pitch angle, defined from the vertical
        aoai = zeros(4,steps); %angle of attack of each external body component(rad)
        phi = zeros(1,steps); %flight path angle
        aoa = zeros(1,steps); %angle of attack

rho = zeros(1,steps); %Nonconstant air density(slug/ft^3)
P_a = zeros(1,steps); %Nonconstant air pressure(lbf/ft^2)
T_a = zeros(1,steps); %Nonconstant air temperature(F)
[rho(1),P_a(1),T_a(1)] = getAtmosphereAlt(0,h0,C1,C2,Temp1);
CP = zeros(1,steps); %CP calculation is different because of aoai

Fx = zeros(1,steps); %Force(total)
Fy = zeros(1,steps);
    % Wy = Weight(lbf) and m = mass(slug) loaded in from Problem_3_new
    % T = zeros(1,steps); T(1) = Ti; %Thrust(lbf)
    T_data(steps) = 0; % pads rest of thrust array with zeros
    Tx = zeros(1,steps);
    Ty = zeros(1,steps); Ty(1) = T_data(1);
    T = zeros(1,steps); T(1) = T_data(1);
    
    Dx = zeros(1,steps); %Drag axial force(lbf)
    Dy = zeros(1,steps);
    CD = zeros(1,steps);
    
    Nx = zeros(1,steps); %Normal force(lbf)
    Ny = zeros(1,steps);
        Nmag = zeros(1,steps); %Magnitude of normal force vector(lbf)
    CNai = zeros(4,steps);
    CN = zeros(1,steps);
if Ty(1)+Wy(1) >= 0 %Don't want negative initial acceleration
    Fy(1) = Ty(1) + Wy(1);
else
    Fy(1) = 0;
end
    
ax = zeros(1,steps); %acceleration(ft/s^2)
ay = zeros(1,steps); ay(1) = Fy(1)/m(1);
a = zeros(1,steps);
    vx = zeros(1,steps); %velocity(ft/s)
    vy = zeros(1,steps);
    v = zeros(1,steps);
        vRWx = zeros(1,steps); %velocity of rocket wrt wind
        vRWy = zeros(1,steps);
        vRW = zeros(1,steps);
        Ma = zeros(1,steps);
        
        vIWx = zeros(4,steps); %velocity of components wrt wind, i elements refer to different components
        vIWy = zeros(4,steps);
        
        vWGx = zeros(1,steps); %velocity of (mean) wind wrt ground
        vWGy = zeros(1,steps);
        vTGx = zeros(1,steps); %velocity of turbulence wrt ground
        vTGy = zeros(1,steps);
    rx = zeros(1,steps); %position(ft)
    rx_max = zeros(1,steps);
    ry = zeros(1,steps);
   
    
%Define arrays at each time slot
for i=2:steps
    
    %Positions and velocities (Euler and such)
    vx(i) = ax(i-1)*dt + vx(i-1);
    vy(i) = ay(i-1)*dt + vy(i-1);
    
    rx(i) = vx(i-1)*dt + rx(i-1);
    rx_max(i) = abs(vx(i-1))*dt + rx_max(i-1);
    ry(i) = vy(i-1)*dt + ry(i-1);
    
    w(i) = aa(i-1)*dt + w(i-1);
    p(i) = w(i-1)*dt + p(i-1);
        
%     %Turbulence
%     % Simulink Variable Initialization
%         h = ry(i);
%         V = hypot(vx(i),vy(i));
%         DCM = angle2dcm(0, p(i), 0); %check angles
% 
%         % Run windTurbulence.slx
%         if i >= 2 && i < steps
%             set_param('windTurbulenceS', 'SimulationCommand', 'start', 'SimulationCommand', 'pause');
%             set_param('windTurbulenceS', 'SimulationCommand', 'continue', 'SimulationCommand', 'pause');
%         elseif i == steps
%             set_param('windTurbulenceS', 'SimulationCommand', 'stop');
%         end
% 
%             Turb_x = out.windTurbulenceSimulink.Data(i,1);
%             Turb_z = out.windTurbulenceSimulink.Data(i,3);

        if ry(i) < 1750
            vTGx(i) = 0;%Turb_x;
            vTGy(i) = 0;%-Turb_z;
        elseif ry(i) >= 1750
            vTGx(i) = 0;%Turb_x*sin(p(i)) + Turb_z*cos(p(i));
            vTGy(i) = 0;%Turb_x*cos(p(i)) - Turb_z*sin(p(i));    
        end
    
    % Wind (starts off the rail)
    %[vTGx(i), vTGy(i)] = windTurbulence(i, steps, ry(i), vx(i), vy(i), p(i));
    %Right now windTurbulence doesn't work properly due to local/global
    %workspace issues with Simulink
    if ry(i)>rail || t(i)>tb
        if ry(i) <= 2000*3.28084
            vWGx(i) = -(windProfileLog(meanWind,rail,ry(i)) + vTGx(i));
        elseif ry(i) > 2000*3.28084
            vWGx(i) = -(windProfilePower(meanWind,rail,ry(i)) + vTGx(i));
        end
        vWGy(i) = vTGy(i);
    end
    
    %Gust at maxQ (tb)
%     if i == 1402
%         gust = 0;%-50; %ft/s
%         vWGx(i) = vWGx(i) + gust;
%     end
    
    %Important derived quantities for later
    v(i) = norm([vx(i), vy(i)]);
    vRWx(i) = vx(i) - vWGx(i);
    vRWy(i) = vy(i) - vWGy(i);
    vRW(i) = norm([vRWx(i), vRWy(i)]);
    
    %Phi calculation
    if vRWy(i)>0
        phi(i) = atan(-vRWx(i)/vRWy(i));
    elseif vRWy(i)==0
        phi(i) = -sign(vRWx(i))*pi/2;
    elseif vRWy(i)<0
        phi(i) = atan(vRWy(i)/vRWx(i)) - sign(vRWx(i))*pi/2;
    else
        error('"You fucked up" -Calvin')
    end
    
    if vy(i)<=0 && t(i)>tb
        p(i) = phi(i);
    end
    aoa(i) = p(i)-phi(i); %opposite sign to actual angle to resolve restoring force sign?
    
%     vIWmag = zeros(4,1);
    for j = 1:4   
        vIWx(j,i) = (CGi(j)-CG(i))*w(i)*cos(p(i)) + vRWx(i);
        vIWy(j,i) = (CGi(j)-CG(i))*w(i)*sin(p(i)) + vRWy(i);
%         vIWmag(j) = norm([vIWx(j,i),vIWy(j,i)]);
        if vIWy(j,i)>0
            phi_prime = atan(-vIWx(j,i)/vIWy(j,i));
        elseif vIWy(j,i)==0
            phi_prime = -sign(vIWx(j,i))*pi/2;
        elseif vIWy(j,i)<0
            phi_prime = atan(vIWy(j,i)/vIWx(j,i)) - sign(vIWx(j,i))*pi/2;
        else
            error('Check')
        end
        aoai(j,i) = p(i) - phi_prime;
    end
    
    %Forces (Weight is done in Problem 3)
    [rho(i),P_a(i),T_a(i)] = getAtmosphereAlt(ry(i),h0,C1,C2,Temp1); %start at h0(FAR altitude)
    
    if t(i)<=tb
        T(i) = T_data(i) + (P_a(1)-P_a(i))*A_e;
        Tx(i) = T(i)*-sin(p(i));
        Ty(i) = T(i)*cos(p(i));
    end
    
    [CD(i), Ma(i)] = get_Cd(t(i),tb,vRW(i),T_a(i),RAS);
    dragProduct = A*CD(i);
%     if vy(i)<=0 && t(i)>tb
%         if ry(i)>mainAlt 
%             dragProduct = drogueArea*drogueCD; %rocket body drag negligible (current model)
%         elseif ry(i)<=mainAlt
%             dragProduct = mainArea*mainCD;
%         end
%     end
    Dx(i) = (0.5*rho(i)*norm([vRWx(i),vRWy(i)])^2)*dragProduct*sin(p(i));
    Dy(i) = (0.5*rho(i)*norm([vRWx(i),vRWy(i)])^2)*dragProduct*-cos(p(i));
        
    CNai0(2) = 4*k*AP*abs(aoai(2,i))/(pi*noseDiameter^2); %Body term
    %dCNaNoseBodyTail = 1; %Safety margin accounting for slender body approximation undershooting
    [CNai(1,i),CNai(2,i),CNai(3,i)] = getCNaNoseBodyTail(Ma(i),noseDiameter,noseLength,bodyDiameter,bodyLength,tailDiameter,tailLength,CNai0,aoai(:,i));
    %CNai(1:3,i) = CNai(1:3,i) + dCNaNoseBodyTail;
    CNai(4,i) = getCNaFins(Ma(i),Ma_dd,finRC,finTC,finSS,finLESA,finsThickness,aoai(4,i),CNai0(4),bodyDiameter);
    CN(i) = sum(CNai(:,i).*aoai(:,i));
    if vy(i)<=0 && t(i)>tb
        CN(i) = 0; %Jank way to disregard lift post-apogee
    end
    
    Nx(i) = (0.5*rho(i)*norm([vRWx(i),vRWy(i)])^2)*A*CN(i)*-cos(p(i));
    Ny(i) = (0.5*rho(i)*norm([vRWx(i),vRWy(i)])^2)*A*CN(i)*-sin(p(i));
    Nmag(i) = norm([Nx(i),Ny(i)]);
    
    %CP
    CP(i) = sum(CNai(:,i).*CPi)/sum(CNai(:,i));

    %Newton's Laws
    aa(i) = sign(Nx(i))*((CP(i)-CG(i))*Nmag(i))/MOI(i);
    ax(i) = (Tx(i)+Nx(i)+Dx(i))/m(i);
    ay(i) = (Ty(i)+Wy(i)+Ny(i)+Dy(i))/m(i);
    if ry(i)==0 && ay(i)<0
        ay(i) = 0;
    end
    a(i) = norm([ax(i), ay(i)]);
    
    %account for ground
    if ry(i)<=0 && t(i)>tb
        tg = t(i);
        tgi = i;
        break
    end
    
    % account for apogee
    if vy(i)<=0 && t(i)>tb
        break
    end
    
    % Stop for emergencies (How to account for nearing apogee???)
%     if abs(aoa(i))>pi/9 && t(i)<tb
%         fprintf('Oh no!\n');
%         fprintf('Time of failure: %.2f\n',(i-1)*dt);
%         fprintf('Value of i at failure: %i\n',i);
%         break
%     end
end
% % End turbulence simulation
% set_param('windTurbulenceS', 'SimulationCommand', 'stop');


%Results

[zmax, tai] = max(ry);
ta = t(tai);
vDeploy_arr(trialNo) = vRW(tai);

rhomax = max(rho);
zb = ry(tb/dt + 1);
vmax = max(sqrt((vx).^2 + (vy).^2));
vb = norm([vx(tb/dt + 1),vy(tb/dt + 1)]);
amax = max(sqrt((ax).^2 + (ay).^2));

for i = 1:steps
    if ry(i)>=rail
        tor = t(i); %time to off-rail
        tori = i;
        break
    end
end

%ORaoa(windex) = aoa(tori);
q = 0.5*rho.*vRW.^2;
[qmax,qmaxi] = max(q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dist_arr(trialNo) = abs(rx(tgi)/5280);
%distT_arr(trialNo) = rx_max(tgi)/5280;
rx_arr(trialNo,:) = rx;
ry_arr(trialNo,:) = ry;
tai_arr(trialNo) = tai;
vMax_arr(trialNo) = vmax;
MaMax_arr(trialNo) = max(Ma);
aMax_arr(trialNo) = max(a);
vOR_arr(trialNo) = v(tori);
aoaOR_arr(trialNo) = aoa(tori);
aoa_arr(trialNo,:) = aoa;
phi_arr(trialNo,:) = phi;
end

%save('FlightSim_apogee_Tsub30','timeLaunch_arr','turbProb_arr','Temp0_arr','meanWind_arr','P0_arr','apogee_arr','vDeploy_arr','aoaOR_arr','aoa_arr');
%save('FlightSim_apogee','timeLaunch_arr','turbProb_arr','Temp0_arr','meanWind_arr','P0_arr','dist_arr','distT_arr','apogee_arr','vDeploy_arr','aoaOR_arr','aoa_arr');
%save('FlightSimWindEffects_T','meanWind_arr','aoaOR_arr','vDeploy_arr','phi_arr');
%save('FlightSimApogeeSpread','tai_arr','rx_arr','ry_arr')

save('FlightSim_dist','dist_arr')

% hold on;
% plot1 = plot(groundWind*0.681818,distMax,'LineWidth',3);
% plot2 = plot(groundWind*0.681818,distMin,'LineWidth',3); 
% 
% x2 = [groundWind*0.681818, fliplr(groundWind*0.681818)];
% inBetween = [distMin, fliplr(distMax)];
% plot3 = fill(x2, inBetween, [0.6 0.6 0.8]);
% 
% xlabel('Windspeed (mph)');
% xlim([0 20]);
% ylabel('Distance (mi)');
% title({'Tranquility 2: Landing Distance Envelope for Various Windspeeds'; '(Total Recovery Failure, Constant Wind Magnitude)'})
% legend([plot1 plot2],{'Maximum Distance','Minimum Distance'});
% legend('Location','northwest');
% grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end
% rx = rx(1:tai);
% 
% ry = ry(1:tai);
% plot(rx,ry);
% xlabel('Distance (ft)(10^3)');
% ylabel('Altitude (ft)(10^3)');
% title('Tranquility 2: Max Altitude = 40,791 ft')
% grid
% axis equal