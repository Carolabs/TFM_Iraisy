function Jacdq = evalJacobiandq(q,qd,SYS)

% Function that evaluates the product of the derivative of the Jacobian 
% matrix of the constraints w.r.t. times the generalized velocities of the
% system
%   q:      Generalized coordinates of the system
%   q:      Generalized velocities of the system
%   SYS:    Structure with system information

x1d = qd(1);
y1d = qd(2);
x2d = qd(3);
y2d = qd(4);
x3d = qd(5);
y3d = qd(6);
sd  = qd(7);

Jacdq = [   2.0*(x1d^2+y1d^2);
            0.0;
            0.0;
            2.0*( (x3d-x2d)^2 + (y3d-y2d)^2 );
            2.0*(x1d^2+y1d^2-sd^2)];
