function HYD = hydSetup()

% This function writes the hydraulic parameters of the system

% _______________________________________________________ System definition
HYD.A   = 0.0065;   % Piston area
HYD.Lc  = 0.442;    % Maximum piston length
HYD.c   = 1.0e5;    % Viscous friction coefficient in actuator
HYD.cd  = 0.67;     % Valve discharge coefficient
HYD.rho = 850.0;    % Fluid density
HYD.pp  = 7.6e6;    % Pump pressure    
HYD.pt  = 0.1e6;    % Tank pressure

% Data from mechanical system
HYD.feff0 = 0.0;      % Effective force from mechanical system
                      % For initialization only

% Bulk modulus parameters
HYD.beta_a = 6.53e-10;     % Taken from Cardona1990
HYD.beta_b = -1.19e-18;

end