% This script simulates the motion of a 2 d.o.f. hydraulically actuated 
% crane using a co-simulation approach.

% Co-simulation scheme:     Non-iterative Jacobi
% Coupling variables:       Force-displacement

function [STORE_HYD, STORE_MECH] = main_Jacobi_fs(hm, hh, st, tfin, corr, coeffs, lim)

addpath('./hyd')
addpath('./mech')


% ___________________________________________________ Simulation parameters

H           = hm;        % Macro step-size
tEnd        = tfin;         % Final simulation time
reportEvery = 5.0e-1;       % Display time in screen to show sim. progress
F_t         = 0;            % Current force
F_t1        = 0;            % Previous force

% ______________________________________________ Create co-simulation units

% First argument: name of subsystem
% Second argument: is hydraulics outputting force? 1 if true, 0 if false
% Third argument (hydraulics only): input for kappa: 0 steps, 1 sinus
hyd     = hydraulics('hyd', 1, 1);
mech    = doublePendulum('mech', 1);

dth = hh;				% Integration step-size in hydraulics
dtm = H;					% Integration step-size in mechanics

hyd.setStepSize(dth);
hyd.setEndTime(tEnd);
hyd.setStoreEvery(st);

mech.setStepSize(dtm);
mech.setEndTime(tEnd);
mech.setStoreEvery(1);

% _________________________________________________________________ Manager
% Processes and stores information

% Select properties and pass MGPROPS as argument during construction
MGPROPS.correctEnergy   = 0;        % Energy correction method
MGPROPS.kp              = 0.5;      % Fraction of dE to be corrected
MGPROPS.kint            = 0.0;      % Fraction of accumulated dE to be corrected
MGPROPS.fCapRel         = 0.05;     % Maximum fraction of fh that can be corrected
MGPROPS.fCapAbs         = 20e3;     % Maximum actuator force cannot exceed this
MGPROPS.dEThold         = H*1;    % Do not correct energy if error < threshold
MGPROPS.xdThold         = 1.0e-3;   % Velocity considered null if below this
MGR = createManager(H,MGPROPS);

% __________________________________________________________ Initialization

% The initialization of the hydraulics requires to know the following
% values that must be provided by the mechanical sub-system
%   - Actuator length s
%   - Actuator length rate sd
%   - Effective force (to evaluate the initial static equilibrium)
% For these reasons, the mechanical sub-system is evaluated first and its
% outputs sent as inputs to the hydraulics sub-system

% Retrieve outputs from the mechanics and send to hydraulics
yMech = mech.readOutputs();
hyd.setInputs(yMech);

% Evaluate effective force of mechanics (for initialization only)
feff0 = mech.evalEffectiveForce;
hyd.setEffectiveForce(feff0);

% Initialize hydraulics
hyd.initialize(0.0);

% Retrieve hydraulics outputs and send to mechanics
yHyd = hyd.readOutputs();
mech.setInputs(yHyd);

% Initialize mechanics
mech.initialize(0.0);

% Initialize manager
yMech = mech.readOutputs;
yHyd = hyd.readOutputs;
MGR = evalManager(MGR, yMech, yHyd, 0.0);

% Correction
F_t = yHyd;
F_t1 = yHyd;

% _______________________________________________________________ Execution
t                   = 0.0;
reportEverySteps    = round(reportEvery/H);
repI                = reportEverySteps - 1;

% The simulation loop follows a Jacobi non-iterative scheme
tic
while t < tEnd

    % Report when required
    repI = repI + 1;
    if (repI == reportEverySteps)
        Xt = ['Macro time: ', num2str(t), '.']; disp(Xt);
        repI = 0;
    end
    
    % Increase goal time
    t = t + H;

    % Exchange coupling variables
    yMech   = mech.readOutputs;
    yHyd    = hyd.readOutputs;
    
    MGR = evalManager(MGR, yMech, yHyd, t);

    % Correction
    F_t1 = F_t;
    F_t = -MGR.uEx(1);

    if corr

        val = coeffs' * [1 F_t1 F_t]';

        if val > lim - MGR.uEx(1)
            val = lim - MGR.uEx(1);
        end

        if val < -lim - MGR.uEx(1)
            val = -lim - MGR.uEx(1);
        end

        mech.setInputs(val);

    else
        mech.setInputs(-MGR.uEx(1));
    end

    
    
    hyd.setInputs(MGR.uEx(2:3));

    % Advance integrations
    doStep(hyd, t);
    doStep(mech, t);

end
toc

% ____________________________________________________________ Post-process
STORE_HYD = hyd.getSTORE();
STORE_MECH = mech.getSTORE();


rmpath('./hyd')
rmpath('./mech')

end
