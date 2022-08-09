function beta = evalbeta(p, SYS)

% Evaluates the bulk modulus of the fluid as a function of beta
%   p:      fluid pressure
%   SYS:    structure with system information

a       = SYS.beta_a;
b       = SYS.beta_b;
beta    = (1.0 + a*p + b*p^2) / (a + 2.0*b*p);

% If beta is taken constant, these reference values can be used 
%beta = 6.1e8;      % Taken from Eryilmaz2006
%beta = 1.15e8;     % Adjusted by trial and error

