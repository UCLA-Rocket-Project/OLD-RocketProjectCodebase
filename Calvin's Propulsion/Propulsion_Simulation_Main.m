clear all; clc;

addpath(genpath('Classes'))
addpath(genpath('Supporting Functions'))
addpath(genpath('Combustion Data'))
load('95EthLOxCompact')

%% Combustion Efficiency
comb_eff = 95;

%% Load Fluid Properties
fluids = Fluids();
gas = fluids.He;

%% Load System Mapping
map = DataMatchAresSystemMapping(fluids);
volumes = VectorizeMapVolumes(map);

%% Simulation Parameters
tf = 20;
dt = 0.001;
n_sim = ceil(tf/dt);
t = zeros(1,n_sim);

%% Array Preparation
ZeroDataArrays(map,n_sim);
imp = zeros(1,n_sim);
i = 0;
while ~OutOfPropellantListener(volumes)
    
    fprintf('i = %d\n',i)
    i = i + 1;
    
    fprintf('%.2f\n',map{1,1}.p)
    fprintf('%.2f\n',map{3,1}.p)
    fprintf('%.2f\n',map{5,1}.p)
    fprintf('%.2f\n',map{5,2}.p)
    fprintf('%.2f\n\n',map{7,1}.p)
    
    fprintf('%.3f\n',map{5,1}.m_prop)
    fprintf('%.3f\n\n',map{5,2}.m_prop)
    
    fprintf('%.3f\n',map{7,1}.FT)
    
    RecordStates(i,map);
    
    Xvals = VectorizeXvals(volumes);
    
    DefineGoverningEquations(map,gas);
    
    A = BlockAMatrix(volumes);

    mdots = MassFlowRates(map,Xvals,A,dt,gas,comb_eff,CombData);
    bvals = mdots;
    
    dX = A\bvals;
    
    Xvals = Xvals + dX*dt;
    
    DealXvals(volumes,Xvals);
    
    imp(i+1) = imp(i) + map{7,1}.data.FT(i)*dt;
    fprintf('%.3f\n\n',imp(i))
    if imp(i+1) >= 8600
        break
    end
    
    t(i+1) = t(i) + dt;
end

%Isp = imp(i)/(map{5,1}.data.m_prop(1) + map{5,2}.data.m_prop(1) - map{5,1}.data.m_prop(i) - map{5,2}.data.m_prop(i))/32.2;

%fprintf('Isp: %.3f\n',Isp)

%% Plotting
lw = 2;
fs = 18;

figure
hold on
plot(t(1:i),map{1,1}.data.p(1:i),'linewidth',lw)
plot(t(1:i),map{3,1}.data.p(1:i),'linewidth',lw)
plot(t(1:i),map{5,1}.data.p(1:i),'linewidth',lw)
plot(t(1:i),map{5,2}.data.p(1:i),'linewidth',lw)
plot(t(1:i),map{7,1}.data.p(1:i),'linewidth',lw)
hold off
legend('Pressurant Tank','Interstitial Tubing','Oxidizer Tank','Fuel Tank','Combustion Chamber')
xlabel('Time (s)','fontsize',fs)
ylabel('Pressure (psia)','fontsize',fs)

figure
plot(t(1:i),map{7,1}.data.FT(1:i),'linewidth',lw)
xlabel('Time (s)','fontsize',fs)
ylabel('Thrust (lbf)','fontsize',fs)

savename = input('file name to save data? ''no'' to not save.\n');

if ~strcmp(savename,'no')
    save(sprintf('%s',savename),'i','t','map','imp','fluids');
end
