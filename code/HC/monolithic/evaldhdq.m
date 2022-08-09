function dhdq = evaldhdq(q,p,qd,kappa,SYS)

% Function that evaluates the partial derivative of the pressure rate
% w.r.t. system coordinates
%   q:      Generalized coordinates of the system
%   p:      Hydraulic pressures
%   qd:     Generalized velocities of the system
%   kappa:  Spool displacement
%   SYS:    Structure with system information

% ______________________________________ Preallocate and retrieve variables
dhdq = zeros(2, 7);

l_1     = 0.5*SYS.Lc + SYS.s0 - q(7);
l_2     = 0.5*SYS.Lc + q(7) - SYS.s0;

% _____________________________________________________ Eval pressure rates
[h1,h2] = evalPressureRates(q,p,qd,kappa,SYS);

% __________________________________________________________________ Update
dhdq(1,7) = h1/l_1;
dhdq(2,7) = -h2/l_2;