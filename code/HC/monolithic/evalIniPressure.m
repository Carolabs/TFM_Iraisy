function [p, kappa, res] = evalIniPressure(p, kappa, q, SYS, niters)

% Function that evaluates the initial pressure p and spool displacement 
% kappa in the system
%
% This requires the solution of a nonlinear system of equations - we solve
% it using a Newton-Raphson approach.
%
%   p:      Hydraulic pressures (initial guess)
%   kappa:  Spool displacement (initial guess)
%   q:      Generalized coordinates of the system (known initial position)
%   SYS:    Structure with system information
%   niters: Number of iterations allowed during the evaluation of initial
%           pressures

% __________________________________________________________ First residual

% Residual evaluation
f = initialization_evalResidual(q, p, kappa, SYS);

% __________________________________________________________ Tangent matrix

dfdz = initialization_evalTangent(q, p, kappa, SYS);

% _______________________________________________________________ Increment
delta_z     = -dfdz\f;

% Update variables
p(1)        = p(1) + delta_z(1);
p(2)        = p(2) + delta_z(2);
kappa       = kappa + delta_z(3);

% _________________________________________________________________ Iterate
for i=0:niters
    
    % Residual evaluation
    f = initialization_evalResidual(q, p, kappa, SYS);
    
    % Tangent matrix evaluation
    dfdz = initialization_evalTangent(q, p, kappa, SYS);
    
    % Evaluate increment
    delta_z     = -dfdz\f;
    
    % Update variables
    p(1)        = p(1) + delta_z(1);
    p(2)        = p(2) + delta_z(2);
    kappa       = kappa + delta_z(3);

end

res = initialization_evalResidual(q, p, kappa, SYS);


end

% _____________________________________________________ Auxiliary functions

% Auxiliary function : evaluation of residual
function f = initialization_evalResidual(q, p, kappa, SYS)

% Residual evaluation
    f       = zeros(3,1);

    qd      = zeros(7,1);
    [h1,h2] = evalPressureRates(q,p,qd,kappa,SYS);

    f(1)    = (2.0*SYS.mp + 2.0*SYS.mh + SYS.m)*SYS.g - (p(2)-p(1))*SYS.A;
    f(2)    = h1;
    f(3)    = h2;

end

% Auxiliary function : evaluation of tangent matrix
function dfdz = initialization_evalTangent(q, p, kappa, SYS)

    % Evaluation of intermediate variables
    p1  = p(1);
    p2  = p(2);

    VALp1   = 2*(SYS.pp-p1)/SYS.rho;
    VALp2   = 2*(SYS.pp-p2)/SYS.rho;
    VAL1T   = 2*(p1-SYS.pt)/SYS.rho;
    VAL2T   = 2*(p2-SYS.pt)/SYS.rho;

    deltap1 = 1.0;
    deltap2 = 1.0;
    deltat1 = 1.0;
    deltat2 = 1.0;

    if (VALp1 < 0); deltap1 = 0; end
    if (VALp2 < 0); deltap2 = 0; end
    if (VAL1T < 0); deltat1 = 0; end
    if (VAL2T < 0); deltat2 = 0; end
    
    l_1     = 0.5*SYS.Lc + SYS.s0 - q(7);
    l_2     = 0.5*SYS.Lc + q(7) - SYS.s0;
    
    beta1 = evalbeta(p1,SYS);
    beta2 = evalbeta(p2,SYS);

    % Tangent matrix
    dfdz        = zeros(3,3);
    dfdz(1,1)   = SYS.A;
    dfdz(1,2)   = -SYS.A;
    dfdz(1,3)   = 0;

    dhdp        = evaldhdp(q, p, kappa, SYS);
    dfdz(2,1)   = dhdp(1,1);
    dfdz(2,2)   = dhdp(1,2);
    dfdz(3,1)   = dhdp(2,1);
    dfdz(3,2)   = dhdp(2,2);

    dfdz(2,3)   =  0.0005*beta1*SYS.cd/(SYS.A*l_1)*...
        (sqrt(VALp1)*deltap1 + sqrt(VAL1T)*deltat1);
    dfdz(3,3)   = -0.0005*beta2*SYS.cd/(SYS.A*l_2)*...
        (sqrt(VALp2)*deltap2 + sqrt(VAL2T)*deltat2);

end