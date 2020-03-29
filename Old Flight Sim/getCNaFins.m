function [finsCNa] = getCNaFins(Ma, Ma_dd, finRC, finTC, finSS, finLESA, finThickness, aoa_fins, finsCNa0, bodyDiameter)
% Returns CNa value for any regime(transsonic, supersonic, even PG compressibility correction)
% Does this account for vortex lift?
% CHECK REFERENCE AREAS

% ACCOUNT FOR DIHEDRAL ANGLES IN 6DOF

    % Model A(refer to http://www.rsandt.com/reports.html supersonic barrowman section)
%     if Ma > sqrt((finRC/finSS)^2 + 1)
% 
%     end

%Subsonic
    % Uses Prandtl-Glauert Compressibility Correction:
    % (http://cambridgerocket.sourceforge.net/AerodynamicCoefficients.pdf)
    if Ma < 0.7
        beta = sqrt(1-Ma^2); % PG correction
        finsCNa = finsCNa0/beta;
    
%Transonic
    elseif Ma>=0.7 && Ma <= 1.2 % Under Ma = 0.7 PG correction still ~valid
        %Temporary transonic placeholder
        finsCNa = finsCNa0;
        % Uses Korn equation: WHY NO WORK??????????
        % (https://pdfs.semanticscholar.org/c754/1e97d55e1e021494050e0ffc38f6702e1554.pdf)
%         K = 0.87; % "Airfoil technology factor", can find from CFD analysis? Default 0.87 for conventional (NACA 6series)
%         t = finsThickness; c = finRC;
%         finsCNa = 10*(K - t/c - Ma_dd)/aoa_fins;
%         %Temporary solution
        
%Supersonic
    elseif Ma > 1.2
    % Uses Busemann second order theory:
    % (http://jestec.taylors.edu.my/eureca2013_4_2014/eureca_13_16_27.pdf)
        beta = sqrt(Ma^2 -1);
        gamma = 1.4; %ratio of specific heats air
        C3 = (gamma*Ma^4 +(Ma^2 -2)^2)/(2*(Ma^2 -1)^1.5);
        A_prime = 2*finThickness/(3*finRC); %thickness effect
        S_wing = (finRC+finTC)*finSS;
        A = (2*finSS)^2 /S_wing; %Aspect ratio
        finsCNa = (4/beta)*(1 -(1 -C3*A_prime)/(2*beta*A));
        
        S_ref = pi*(bodyDiameter/2)^2;
        K_wb = 1 + ((bodyDiameter/2)/finSS)^1.15;
        finsCNa = K_wb*(S_wing/S_ref)*finsCNa;
    % Uses flat plate linear theory
    % Uses following source:
    % (http://e.roohi.profcms.um.ac.ir/imagesm/1019/stories/PDFs/Supersonic/pages%20from%20bertin.pdf)
    % (Pasquale Manned Spacecraft Design Principles Eq. 7.138, 7.139)
%         % Important Variables
%         beta = sqrt(Ma^2 - 1);
%         S_fin = (finRC+finTC)*finSS;
%         S_ref = pi*(bodyDiameter/2)^2;
%         K_wb = 1 + ((bodyDiameter/2)/finSS)^1.15;
%         finsCNa = K_wb*(S_fin/S_ref)*4/beta;
        
    % Uses NACA TN 1977:
    % (https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19930082644.pdf)
%         % Important Variables
%         beta = sqrt(Ma^2 - 1);
%         lambda = finTC/finRC;
%         A = (2*finSS)^2 / ((finRC+finTC)*finSS);
%         m = cot(finLESA);
%         w = 4*m/(A*(1+lambda));
%         n = 1 + w*(lambda-1);
%         k = sqrt(1 - (beta*m)^2);
%         syms z; expr = sqrt(1 - (k*sin(z))^2);
%         E = int(expr,[0 pi/2]);
% 
%         if finLESA > pi/2 - asin(1/Ma)
%         %Subsonic leading edge (inside Mach cone)
%         % taper ratio(lambda), beta, cotfinLESA(m), aspect ratio(A),w,n
% 
%             %Limiting Condition
%             test1 = beta*A*(1+lambda);
%             test2 = test1/(test1 +4*(1-lambda));
%             test3 = test1/(4-test1);
%             
%             if test1 >= 2
%                 assert(test2<=beta*m && beta*m<=1);
%             elseif test1 <2
%                 assert(test2<=beta*m && beta*m<=test3);
%             else
%                 error('test1 in getCNaFins effed up');
%             end
%             
%             % Part1 (coefficient of PartA)
%             Part1 = A/E;
%             % PartA (first big parentheses) has 3 parts (a,b,c)
%             PartA_a_1 = w^2 / (1-n^2)^3/12; % IS 3/12 A TYPO IN THE PAPER??????
%             PartA_a_2 = asin(n) - asin(((beta*m+1)*(n^2-1)+w*(beta*m*n+1)) / (w*(beta*m+n)));
%             PartA_b_a_1 = (1+beta*m)^0.5 / (1-beta*m)^1.5;
%             PartA_b_a_2 = acos((1 +beta*m*n +w*(beta*m -1))/(beta*m +n));
%             PartA_b_b = n*w^2 /(1 -n^2);
%             PartA_c_1 = (w*n*(beta*m -1) +beta*m*(n^2 -1))*(1 +beta*m)^0.5 / ((beta*m +n)*(n^2 -1)*(beta*m -1));
%             PartA_c_2 = ((w+n-1)*((beta*m +1)*(n+1) +w*(beta*m -1)))^0.5;
%             % Part2 (coefficient of PartB)
%             Part2 = 4*A/(pi*(1 +beta*m)^0.5);
%             % PartB (second big parentheses) has 3 parts (a,b,c)
%             PartB_a_1 = (1+n+w)^2 /(4*(1+n)^1.5);
%             PartB_a_2 = acos(((n+w)*(beta*m -n) +2*(1-w) +beta*m +n) / ((1+n+w)*(beta*m +n)));
%             PartB_b_1 = -1 / (1 - beta*m)^1.5;
%             PartB_b_2 = acos((1 +beta*m*n +w*(beta*m -1)) / (beta*m +n));
%             PartB_c_1 = ((1 +beta*m)*(1+n) -w*(1 -beta*m)) / (2*(beta*m +n)*(1 -beta*m)*(1+n));
%             PartB_c_2 = ((w-n+1)*((beta*m +1)*(n+1) +w*(beta*m -1)))^0.5;
%         
%             PartA_a = PartA_a_1*PartA_a_2;
%             PartA_b = (PartA_b_a_1*PartA_b_a_2) +PartA_b_b;
%             PartA_c = PartA_c_1*PartA_c_2;
%             PartB_a = PartB_a_1*PartB_a_2;
%             PartB_b = PartB_b_1*PartB_b_2;
%             PartB_c = PartB_c_1*PartB_c_2;
%             
%             PartA = PartA_a+PartA_b+PartA_c;
%             PartB = PartB_a+PartB_b+PartB_c;
%             
%             finCNa = Part1*PartA + Part2*PartB;
%             
%         elseif finLESA <= pi/2 - asin(1/Ma)
%         %Supersonic leading edge (outside Mach cone)
%         % beta,cotfinLESA(m),aspect ratio(A)
%         
%             %Limiting Condition
%             assert(beta > 2/A + (1+n)/(2*m));
%         
%             
%             Part1 = 1 / ((beta*m)^2 -1)^0.5;
%             Part2_a_1 = (4*beta*m -beta*A*n +beta*A)^2 /(2*pi*beta*A*(1-n^2));
%             Part2_a_2_a = n*acos(1/(beta*m));
%             Part2_a_2_b = ((beta*m)^2 -1)^0.5 * acos(-n/(beta*m)) / ((beta*m)^2 - n^2);
%             Part2_b = -(4*beta*m +beta*A*n -beta*A)^2 * (beta*m +1)^0.5 / (4*beta*A*(1-n)*(beta*m +n)^0.5);
%             
%             finCNa = (1/beta)*Part1*(Part2_a_1*(Part2_a_2_a+Part2_a_2_b) + Part2_b);
%   
%         end
    else
        error('Mach Number messed up: getCNaFins');
    end