function q = ValveChokedGasFlow(Cv,p_1,T_1,G_g)

N2 = 22.67;
q = 0.471*N2*Cv*p_1*sqrt(1/(G_g*T_1));