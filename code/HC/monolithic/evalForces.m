function Q = evalForces(p, qd, SYS)

% Function that evaluates the generalized applied forces
%   p:      Hydraulic pressures
%   qd:     Generalized velocities of the system
%   SYS:    Structure with system information

% ______________________________________________________ Retrieve variables

m   = SYS.m;
g   = SYS.g;
mp  = SYS.mp;
mh  = SYS.mh;
A   = SYS.A;
c   = SYS.c;

sd  = qd(7);

p1  = p(1);
p2  = p(2);

% __________________________________________________________________ Update
Q = [0.0; -m*g; 0.0; -mp*g; 0.0; -mh*g; (p2-p1)*A - c*sd];