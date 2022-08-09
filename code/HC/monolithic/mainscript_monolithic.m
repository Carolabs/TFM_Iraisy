% This script performs the numerical integration of a double pendulum with
% hydraulic actuation. The approach is similar to the one introduced in:
%
% M. A. Naya, J. Cuadrado, D. Dopico, U. Lugrís
% An efficient unified method for the combined simulation of multibody and
% hydraulic dynamics: comparison with simplified and co-integration
% approaches
% The Archive of Mechanical Engineering, 58(2):223-243, 2011.
% d.o.i.: 10.2478/v10180-011-0016-4
%
% The system dynamics is solved following a monolithic approach (mechanics
% and hydraulics are solved together)
%
% The pendulum is composed of a rod with a uniformly distributed mass and a
% massless rod, plus two point masses, one at the end of each rod.
%
% The mechanical system is modelled with the x and y coordinates of points 
% 1 (centre of first rod), 2 (tip of first rod), and 3 (tip of second rod).
% The actuator length s is also used as coordinate.
% The mass matrix is singular, but this does not seem to bother the
% algorithm too much.

function STORAGE = mainscript_monolithic(tend, step)

% __________________________________________________________ Output options
OPT.showAnimation   = 0;
OPT.showPlots       = 1;
OPT.saveToFile      = 1;

% ___________________________________________________ Simulation properties
saveRate    = 1;

% Method parameters
SYS.dt      = step;       % Integration step-size
SYS.maxIter = 100;          % Maximum number of Newton-Raphson iteration
SYS.maxError= 1.0e-5;       % Maximum error in Newton-Raphson iteration
SYS.nProjs  = 3;            % Number of velocity and accel. projections
SYS.tEnd    = tend;         % Final simulation time
SYS.alpha   = 1.0e12;       % Penalty factor for augmented Lagrangian meth.
SYS.xi      = 1.0;          % Stabilization parameters (initial acc.)
SYS.w       = 10.0;

SYS.inputK  = 1;            % 0: steps; 1: sin

% ______________________________________________________________ Properties

SYS.L   = 1.0;          % First rod length
SYS.m   = 200.0;        % First rod mass (distributed)
SYS.mp  = 250.0;        % Point mass at the end of first rod
SYS.mh  = 100.0;        % Point mass at the end of second rod
SYS.g   = 9.81;         % Gravity
SYS.L23 = SYS.L/2.0;    % Length of second rod

SYS.xA  = 0.0;
SYS.yA  = 0.0;
SYS.xB  = sqrt(3.0)/2.0;    % x coordinate of fixed point B
SYS.yB  = 0.0;
SYS.th  = deg2rad(30.0);   % Initial angle

SYS.A   = 0.0065;   % Piston area
SYS.Lc  = 0.442;    % Maximum piston length
SYS.c   = 1.0e5;    % Viscous friction coefficient in actuator
SYS.cd  = 0.67;     % Valve discharge coefficient
SYS.rho = 850.0;    % Fluid density
SYS.pp  = 7.6e6;    % Pump pressure    
SYS.pt  = 0.1e6;    % Tank pressure
SYS.s0  = 0.5;      % Initial length of actuator

% Bulk modulus parameters
SYS.beta_a = 6.53e-10;     % Taken from Cardona1990
SYS.beta_b = -1.19e-18;


% _________________________________________________________________ Storage

% Number of storage points
npoints     = round(max(SYS.tEnd/(SYS.dt*saveRate), 1));

% Store time
nvars           = 7;
nconst          = 5;
STORAGE.t       = zeros(1, npoints);
STORAGE.kappa   = zeros(1, npoints);
STORAGE.pos     = zeros(nvars, npoints);
STORAGE.vel     = zeros(nvars, npoints);
STORAGE.acc     = zeros(nvars, npoints);
STORAGE.p       = zeros(2, npoints);
STORAGE.pd      = zeros(2, npoints);
STORAGE.lambda  = zeros(nvars - 2, npoints);
STORAGE.F       = zeros(1, npoints);
STORAGE.violC   = zeros(1, npoints);

% ___________________________________________________ Initial configuration
L  = SYS.L;
th = SYS.th;
s0 = SYS.s0;
q  = [  L/2.0*cos(th); L/2.0*sin(th); ...
        L*cos(th); L*sin(th); ...
        L*cos(th); L*sin(th) - SYS.L23; s0];

% Initial velocities are set to zero
qd  = zeros(nvars,1);

% Evaluate initial pressures
p       = [3.3e6;4.4e6];        % Initial guess
kappa   = 0.5;                  % Spool displacement (valve)
[p, kappa, res] = evalIniPressure(p, kappa, q, SYS, 4);

% Compute Initial Force
F = (p(2) - p(1))*SYS.A;

% Store initial value of the spool displacement
SYS.kappa0 = kappa;


% % ________________________________________________ Dynamic terms (constant)

% Mass matrix
MM  = [ 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0;
        0, 0, SYS.mp+SYS.m/3, 0, 0, 0, 0;
        0, 0, 0, SYS.mp+SYS.m/3, 0, 0, 0;
        0, 0, 0, 0, SYS.mh, 0, 0;
        0, 0, 0, 0, 0, SYS.mh, 0;
        0, 0, 0, 0, 0, 0, 0];
    
% Damping matrix
C               = zeros(nvars,nvars);
C(nvars,nvars)  = SYS.c;

% Partial derivative of forces w.r.t. pressures
dQdp        = zeros(nvars,2);
dQdp(7,1)   = -SYS.A;
dQdp(7,2)   = SYS.A;

% ________________________________________________ Dynamic terms (variable)

% Kinematic constraints
Phi     = evalConstr(q, SYS);
Jac     = evalJacobian(q, SYS);
Jacdq   = evalJacobiandq(q, qd, SYS);
 
% Partial derivative terms
dhdq    = evaldhdq(q,p,qd,kappa,SYS);
dhdqd   = evaldhdqd(q,p,SYS);

% ___________________________________________________ Initial accelerations

% Evaluate pressure rates (should be zero)
[h1,h2] = evalPressureRates(q,p,qd,kappa,SYS);
pd = [h1;h2];

% Applied forces
Q   = evalForces(p, qd, SYS);

% Evaluate accelerations
LEAD_iniacc = [MM Jac'; Jac, zeros(5,5)];
RHS_iniacc  = [Q; -Jacdq];
sol_iniacc  = LEAD_iniacc\RHS_iniacc;
qdd         = sol_iniacc(1:7);
lambda      = sol_iniacc(8:12);

% ______________________________________________________ Dynamic simulation

% Initialize
t           = 0.0;              % Time
i           = saveRate-1;       % Storage counters
storeIdx    = 0;


while (t<=SYS.tEnd)
    
  	% Save values
    i = i+1;
    if (i == saveRate)
        
        t
        storeIdx = storeIdx + 1;
        
        STORAGE.t(storeIdx)         = t;
        STORAGE.kappa(storeIdx)     = kappa;
        STORAGE.pos(:,storeIdx)     = q;  
        STORAGE.vel(:,storeIdx)     = qd; 
        STORAGE.acc(:,storeIdx)     = qdd;
        STORAGE.lambda(:,storeIdx)  = lambda;
        STORAGE.F(storeIdx)         = F;
        STORAGE.violC(storeIdx)     = norm(Phi,2);
        STORAGE.p(:,storeIdx)       = p;
        STORAGE.pd(:,storeIdx)      = pd;
        
        i = 0;
    end

    % Predictor
    qdOldBow    = -( 2.0/SYS.dt * q +  qd );
    qddOldBow   = -( 4.0/SYS.dt^2 * q  + 4.0/SYS.dt * qd + qdd);
    pdOldBow    = -( 2.0/SYS.dt * p +  pd ); 

    q_n1        = q + SYS.dt * qd + (SYS.dt^2/2.0) * qdd;
    qd_n1       = qd + SYS.dt * qdd;
    qdd_n1      = qdd;
    p_n1        = p + SYS.dt * pd;
    pd_n1       = pd;
    lambda_n1   = lambda;

    % Corrector
    n_iter = 0;
    while n_iter <= SYS.maxIter

        % Re-evaluate dynamic terms
        Phi     = evalConstr(q_n1, SYS);
        Jac     = evalJacobian(q_n1, SYS);
        Jacdq   = evalJacobiandq(q_n1, qd_n1, SYS);
        Q       = evalForces(p_n1, qd_n1, SYS);
        [h1,h2] = evalPressureRates(q_n1,p_n1,qd_n1,kappa,SYS);
        h       = [h1;h2];
        dhdq    = evaldhdq(q_n1,p_n1,qd_n1,kappa,SYS);
        dhdqd   = evaldhdqd(q_n1,p_n1,SYS);
        dhdp    = evaldhdp(q_n1, p_n1, kappa, SYS);
        
        % Update Lagrange multipliers
        lambda_n1 = lambda_n1 + SYS.alpha*Phi; 

        % Eval residual
        f = SYS.dt^2*0.25*...
            [MM*qdd_n1 + Jac'*SYS.alpha*Phi + Jac'*lambda_n1 - Q;
             pd_n1 - h];

        % Eval tangent matrix
        dfdx11 = MM + SYS.dt/2.0*C + SYS.dt^2/4.0*(Jac'*SYS.alpha*Jac);
        dfdx12 = -SYS.dt^2/4.0 * dQdp;
        dfdx21 = -SYS.dt/2.0*(SYS.dt/2.0*dhdq + dhdqd);
        dfdx22 = SYS.dt/2.0 * (eye(2) - SYS.dt/2.0*dhdp);
        dfdx    = [dfdx11, dfdx12; dfdx21, dfdx22];

        % Evaluate increment
        dx = -dfdx\f;

        % Update positions, velocities and accelerations, and pressures
        q_n1    = q_n1 + dx(1:7);
        qd_n1   = qdOldBow + 2.0/SYS.dt * q_n1;
        qdd_n1  = qddOldBow + 4.0/SYS.dt^2 * q_n1;

        p_n1    = p_n1 + dx(8:9);
        pd_n1   = pdOldBow + 2.0/SYS.dt * p_n1; 

        % Check error
        error_e = norm(dx, 2);
        if abs(error_e) < abs(SYS.maxError)
            break;     
        end 

        % Increase number of iterations, warn if reached max iterations
        n_iter = n_iter+1;
        if (n_iter == SYS.maxIter-1)
            disp ('Reached max. num. of iterations');
        end
    end
    
    % Increase time and move to next integration time step
    t       = t + SYS.dt;
    q       = q_n1;
    qd      = qd_n1;
    qdd     = qdd_n1;
    p       = p_n1;
    F       = (p(2) - p(1))*SYS.A - SYS.c*qd(7);
    pd      = pd_n1;
    lambda  = lambda_n1;
    
    % Need to project velocities and accelerations
    W           = MM + SYS.dt/2.0 * C;
    ProjLead    = W + SYS.dt^2/4.0*Jac'*SYS.alpha*Jac;
    
    for j=1:SYS.nProjs
        ProjRHS     = [W*qd, ...
                        W*qdd - SYS.dt^2/4.0*Jac'*SYS.alpha*Jacdq];
        solProj     = ProjLead\ProjRHS;

        qd          = solProj(:,1);
        qdd         = solProj(:,2);
    end
    
    % Update the valve parameter
    kappa       = updateKappa(t, SYS);
    
end


end