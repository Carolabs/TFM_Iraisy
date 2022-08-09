function A = evalJacobian(q, SYS)

% Function that evaluates the Jacobian matrix of the constraints
%   q:      Generalized coordinates of the system
%   SYS:    Structure with system information

% ______________________________________ Retrieve variables and preallocate
A = zeros(5,7);

xA = SYS.xA;
yA = SYS.yA;
xB = SYS.xB;
yB= SYS.yB;

x1 = q(1);
y1 = q(2);
x2 = q(3);
y2 = q(4);
x3 = q(5);
y3 = q(6);
s  = q(7);

% ________________________________________________________________ Evaluate
A(1,1) = 2.0*(x1-xA);
A(1,2) = 2.0*(y1-yA);
A(2,1) = -2.0;
A(2,3) = 1.0;
A(3,2) = -2.0;
A(3,4) = 1.0;
A(4,3) = -2.0*(x3-x2);
A(4,4) = -2.0*(y3-y2);
A(4,5) =  2.0*(x3-x2);
A(4,6) =  2.0*(y3-y2);
A(5,1) = 2.0*(x1 - xB);
A(5,2) = 2.0*(y1 - yB);
A(5,7) = -2.0*s;