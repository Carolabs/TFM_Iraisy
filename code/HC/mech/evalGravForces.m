function fg = evalGravForces(th, MECH)

% Evaluates the gravitational forces of the mechanical system

% _________________________________________________________ Retrieve angles
th1 = th(1);
th2 = th(2);

% ________________________________________________________________ Evaluate
fg(1,1) = -cos(th1) * MECH.L * MECH.g * (0.5*MECH.m + MECH.mp + MECH.mh);
fg(2,1) = - MECH.mh * MECH.L23 * MECH.g * cos(th2);