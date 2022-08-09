function [p, kappa, res] = evalIniPressure(p, kappa, s, HYD, feff0, niters)

% Function that evaluates the initial pressure p and spool displacement 
% kappa in the system
%
% This requires the solution of a nonlinear system of equations - we solve
% it using a Newton-Raphson approach.
%
%   p:      Hydraulic pressures (initial guess)
%   kappa:  Spool displacement (initial guess)
%   s:      Cylinder displacement (initially known)
%   HYD:    Structure with hydraulics information
%   feff0:  Initial force exerted by the cylinder
%   niters: Number of iterations allowed during the evaluation of initial
%           pressures

% __________________________________________________________ First residual

% Residual evaluation
f = initialization_evalResidual(s, p, kappa, HYD, feff0);

% __________________________________________________________ Tangent matrix

dfdz = initialization_evalTangent(s, p, kappa, HYD);

% _______________________________________________________________ Increment
delta_z     = -dfdz\f;

% Update variables
p(1)        = p(1) + delta_z(1);
p(2)        = p(2) + delta_z(2);
kappa       = kappa + delta_z(3);

% _________________________________________________________________ Iterate
    for i=0:niters

        % Residual evaluation
        f = initialization_evalResidual(s, p, kappa, HYD, feff0);

        % Tangent matrix evaluation
        dfdz = initialization_evalTangent(s, p, kappa, HYD);

        % Evaluate increment
        delta_z     = -dfdz\f;

        % Update variables
        p(1)        = p(1) + delta_z(1);
        p(2)        = p(2) + delta_z(2);
        kappa       = kappa + delta_z(3);

    end

res = initialization_evalResidual(s, p, kappa, HYD, feff0);


end

% _____________________________________________________ Auxiliary functions

% Auxiliary function : evaluation of residual
function f = initialization_evalResidual(s, p, kappa, HYD, feff0)

% Residual evaluation
    f       = zeros(3,1);

    [h1,h2] = evalPressureRates(s, 0.0, p,kappa,HYD);
    
    % Here we are using the initial effective force provided by the
    % mechanical system
    f(1)    = -feff0 - (p(2)-p(1))*HYD.A;
    f(2)    = h1;
    f(3)    = h2;

end

% Auxiliary function : evaluation of tangent matrix
function dfdz = initialization_evalTangent(s, p, kappa, HYD)

    % Evaluation of intermediate variables
    p1  = p(1);
    p2  = p(2);

    VALp1   = 2*(HYD.pp-p1)/HYD.rho;
    VALp2   = 2*(HYD.pp-p2)/HYD.rho;
    VAL1T   = 2*(p1-HYD.pt)/HYD.rho;
    VAL2T   = 2*(p2-HYD.pt)/HYD.rho;

    deltap1 = 1.0;
    deltap2 = 1.0;
    deltat1 = 1.0;
    deltat2 = 1.0;

    if (VALp1 < 0); deltap1 = 0; end
    if (VALp2 < 0); deltap2 = 0; end
    if (VAL1T < 0); deltat1 = 0; end
    if (VAL2T < 0); deltat2 = 0; end
    
    l_1     = 0.5*HYD.Lc + HYD.s0 - s;
    l_2     = 0.5*HYD.Lc + s - HYD.s0;
    
    beta1 = evalbeta(p1,HYD);
    beta2 = evalbeta(p2,HYD);

    % Tangent matrix
    dfdz        = zeros(3,3);
    dfdz(1,1)   = HYD.A;
    dfdz(1,2)   = -HYD.A;
    dfdz(1,3)   = 0;

    dhdp        = evaldhdp(s, p, kappa, HYD);
    dfdz(2,1)   = dhdp(1,1);
    dfdz(2,2)   = dhdp(1,2);
    dfdz(3,1)   = dhdp(2,1);
    dfdz(3,2)   = dhdp(2,2);

    dfdz(2,3)   =  0.0005*beta1*HYD.cd/(HYD.A*l_1)*...
        (sqrt(VALp1)*deltap1 + sqrt(VAL1T)*deltat1);
    dfdz(3,3)   = -0.0005*beta2*HYD.cd/(HYD.A*l_2)*...
        (sqrt(VALp2)*deltap2 + sqrt(VAL2T)*deltat2);

end