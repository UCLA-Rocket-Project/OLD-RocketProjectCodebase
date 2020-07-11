function q = ValveUnchokedGasFlow(Cv,p_1,p_2,T_1,G_g)

N2 = 22.67;
dp = p_1 - p_2;
q = N2*Cv*p_1*(1 - 2*dp/(3*p_1))*sqrt(dp/(p_1*G_g*T_1));