function [s, sd, Ai] = findSfromAngles(th,thd,MECH)

% Returns the actuator length, its rate, and the Jacobian matrix that
% relates sd and the angular velocities of the mechanism
%   th:     Joint angles
%   thd:    Joint rates
% 	MECH: 	System properties
%   s:      Cylinder length
%   sd:     Cylinder length rate
%   Ai:     Jacobian matrix

% _________________________________________ Retrieve angles and angle rates
th1     = th(1);
th1d    = thd(1);
th2d    = thd(2);

% ________________________________________________________________ Evaluate

% Cylinder length
s = MECH.L * sqrt(1.0 - sqrt(3.0)/2 * cos(th1));

% Jacobian matrix
Ai = [sqrt(3.0)*MECH.L^2/4.0 * sin(th1)/s, 0];

% Cylinder length rate
sd = Ai * [th1d;th2d];