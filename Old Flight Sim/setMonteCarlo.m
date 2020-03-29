function [timeLaunch, turbProb, T0, meanWind, P0] = setMonteCarlo(N_n, mu_all, sigma_all, no_time, no_turbProb)
% Insert description

    %% Set uniform variable ranges
    turbProb = no_turbProb; % (three cases: 10^-1, 10^-2, 10^-3)
    no_timeLaunch = no_time; % (units of 10 min from 7-12 AM)
    
    %% Get initial condition array (run normal Latin hypercube at each turb/time combo)
    N_u = turbProb*no_timeLaunch;
    rows = N_n*N_u; % Total number of random sets of variables
    initCond = zeros(rows, 5);
    i_N_u = 1;
    for i_timeLaunch = 1:no_timeLaunch
        for i_turbProb = 1:turbProb
            initCond(((i_N_u-1)*N_n + 1):i_N_u*N_n,1) = i_timeLaunch;
            initCond(((i_N_u-1)*N_n + 1):i_N_u*N_n,2) = i_turbProb;
            initCond(((i_N_u-1)*N_n + 1):i_N_u*N_n,3:5) = lhsnorm(mu_all(:,i_timeLaunch),sigma_all(:,:,i_timeLaunch),N_n);
            i_N_u = i_N_u + 1;
        end
    end    
     
    %% Setup Output Variables
    timeLaunch = initCond(:,1);
    turbProb_range = initCond(:,2);
    T0 = initCond(:,3);
    meanWind = initCond(:,4);
    P0 = initCond(:,5);
    
    %% Fixing turbProb (to output proper strings for Simulink input)
    turbProb = strings(rows,1);
    for i = 1:rows
        switch turbProb_range(i)
            case 1
                turbProb(i) = "10^-1"; % Very light turbulence at high altitudes
            case 2
                turbProb(i) = "10^-2"; % Light turbulence at high altitudes
            case 3
                turbProb(i) = "10^-3"; % Moderate turbulence at high altitudes
            otherwise
                error('Boundary condition problem!');
        end
    end
    
end