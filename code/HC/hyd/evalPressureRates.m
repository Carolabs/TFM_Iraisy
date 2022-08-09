function [h1,h2] = evalPressureRates(s, sd, p, kappa, HYD)

% Evaluates the pressure rates of the system
%   s:      Cylinder length
%   sd:     Cylinder length rate
%   p:      Hydraulic pressures
%   kappa:  Spool displacement
%   HYD:    Structure with hydraulics information

% ___________________________________________________________ Retrieve vars
A       = HYD.A;
Lc      = HYD.Lc;
cd      = HYD.cd;
pp      = HYD.pp;
pt      = HYD.pt;
rho     = HYD.rho;
s0      = HYD.s0;

p1      = p(1); 
p2      = p(2); 

% __________________________________________________ Intermediate variables

% Evaluate lengths
l_1     = 0.5*Lc + s0 - s;
l_2     = 0.5*Lc + s - s0;

% Evaluate areas
Ai      = 0.0005*kappa;
Ao      = 0.0005*(1.0-kappa);

VALp1   = 2*(pp-p1)/rho;
VALp2   = 2*(pp-p2)/rho;
VAL1T   = 2*(p1-pt)/rho;
VAL2T   = 2*(p2-pt)/rho;

deltap1 = 1.0;
deltap2 = 1.0;
deltat1 = 1.0;
deltat2 = 1.0;

if (VALp1 < 0); deltap1 = 0; end
if (VALp2 < 0); deltap2 = 0; end
if (VAL1T < 0); deltat1 = 0; end
if (VAL2T < 0); deltat2 = 0; end

beta1 = evalbeta(p1,HYD);
beta2 = evalbeta(p2,HYD);


% __________________________________________________________ Pressure rates
h1 = beta1/(A*l_1) * ...
    ( A*sd + Ai*cd*sqrt(VALp1)*deltap1 - Ao*cd*sqrt(VAL1T)*deltat1);
h2 = beta2/(A*l_2) * ...
    (-A*sd + Ao*cd*sqrt(VALp2)*deltap2 - Ai*cd*sqrt(VAL2T)*deltat2);

