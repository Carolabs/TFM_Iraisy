function MECH = mechSetup()

% This function writes the mechanical parameters of the system

% _______________________________________________________ System definition

MECH.L   = 1.0;         % First rod length
MECH.L23 = MECH.L/2.0;  % Second rod length
MECH.m   = 200.0;       % First rod mass (distributed)
MECH.mp  = 250.0;       % Point mass at the end of first rod
MECH.mh  = 100.0;       % Point mass at the end of second rod
MECH.g   = 9.81;        % Gravity

% Coordinates of fixed points
MECH.xA  = 0.0;
MECH.yA  = 0.0;
MECH.xB  = sqrt(3.0)/2.0;    % x coordinate of fixed point B
MECH.yB  = 0.0;

end