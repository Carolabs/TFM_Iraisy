function dhdqd = evaldhdqd(q,p,SYS)

% Function that evaluates the partial derivative of the pressure rate
% w.r.t. system velocities
%   q:      Generalized coordinates of the system
%   p:      Hydraulic pressures
%   SYS:    Structure with system information

% ______________________________________ Preallocate and retrieve variables
dhdqd = zeros(2, 7);

l_1     = 0.5*SYS.Lc + SYS.s0 - q(7);
l_2     = 0.5*SYS.Lc + q(7) - SYS.s0;

p1 = p(1);
p2 = p(2);

% _____________________________________________________ Eval pressure rates
beta1 = evalbeta(p1, SYS);
beta2 = evalbeta(p2, SYS);

% __________________________________________________________________ Update
dhdqd(1,7) = beta1/l_1;
dhdqd(2,7) = -beta2/l_2;