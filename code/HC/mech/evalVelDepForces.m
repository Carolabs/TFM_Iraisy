function fc = evalVelDepForces(th, thd, MECH)

% Evaluates the velocity dependent and Coriolis term of the mechanical system

% _________________________________________ Retrieve angles and angle rates
th1     = th(1);
th2     = th(2);
th1d    = thd(1);
th2d    = thd(2);

% ________________________________________________________________ Evaluate
fc = MECH.L * MECH.L23 * MECH.mh * sin(th1-th2) * [th2d^2; - th1d^2];