function [vTGx,vTGy] = windTurbulence(i, steps, ry, vx, vy, pitch)
% 3DOF SIM: Runs the windTurbulence Simulink file. dt = 0.01s
% windTurbulence.slx currently uses Von Karman Continuous Model (+q -r)
% Note: MUST START AND END SIMULINK MODEL OUTSIDE OF THIS FUNCTION!!!

% Input: 
% Output: 

% Notes: add set dt, wingspan functionality
%% Simulink Variable Initialization
h = ry;
V = hypot(vx,vy);
DCM = angle2dcm(0, pitch, 0); %check angles

%% Run windTurbulence.slx
if i >= 2 && i < steps
    set_param('windTurbulenceS', 'SimulationCommand', 'start', 'SimulationCommand', 'pause');
    set_param('windTurbulenceS', 'SimulationCommand', 'continue', 'SimulationCommand', 'pause');
elseif i == steps
    set_param('windTurbulenceS', 'SimulationCommand', 'stop');
end

    Tx = out.windTurbulenceSimulink.Data(i,1);
    Tz = out.windTurbulenceSimulink.Data(i,3);
    
if ry < 1750
    vTGx = Tx;
    vTGy = -Tz;
elseif ry > 1750
    vTGx = Tx*sin(pitch) + Tz*cos(pitch);
    vTGy = Tx*cos(pitch) - Tz*sin(pitch);    
end
end