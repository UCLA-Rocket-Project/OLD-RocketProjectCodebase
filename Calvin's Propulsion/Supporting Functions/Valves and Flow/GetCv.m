function Cv = GetCv(q,p1,p2,G_g,T1,gam_g)

Cv_guess = 0;
q_guess = 0;

while abs(q_guess - q) > 1E-1
    Cv_guess = Cv_guess + 0.0001;
    q_guess = ValveGasFlow(Cv_guess,p1,p2,G_g,T1,gam_g);
end

Cv = Cv_guess;