function MM = evalMassMatrix(th, MECH)

% Evaluates the mass matrix of the mechanical system

% _________________________________________________________ Retrieve angles
th1 = th(1);
th2 = th(2);

% ____________________________________________________ Terms in mass matrix
mth1    = MECH.L^2 * (1.0/3.0 * MECH.m + MECH.mp + MECH.mh);
mth2    = MECH.L23^2 * MECH.mh;
mcp     = MECH.L * MECH.L23 * cos(th1-th2) * MECH.mh;

MM = [mth1, mcp; mcp, mth2];