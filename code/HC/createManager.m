function MGR = createManager(H, PR)

% Create the manager to handle energy monitoring of hydraulic crane
% H:                Macro step-size
% PR:               Properties
%   .correctEnergy  Energy correction method. 
%                       0: no correction
%                       1: remove dEk/2

% Create structures for storage
MGR.STORE       = [];       
MGR.storeIdx    = 0;
MGR.uTilde1     = 0.0;      % Required for power indicator
MGR.uTilde2     = 0.0;
MGR.deltaP      = 0.0;      % Residual power
MGR.deltaE      = 0.0;      % Residual energy
MGR.Ek          = 0.0;      % Accumulated residual energy

% Corrections
MGR.fCorr       = 0.0;      % Corrective force
MGR.dWfCorr     = 0.0;      % Work of the corrective force during a step
MGR.Ebalance    = 0.0;      % Energy balance: (deltaE/2.0 - dWfCorr)

% _________________________________________________________________________
%                                     Integration and correction properties
MGR.H               = H;
MGR.PR              = PR;

% _________________________________________________________________________
%                                                       Connectivity matrix
MGR.L = [ 0 0 -1; 
          1 0 0;
          0 1 0];

% _________________________________________________________________________
%                                                         Transition matrix
MGR.Z1      = zeros(3);
MGR.Z       = zeros(3);
MGR.rhoInf  = 0;
MGR.detZ    = 0;
MGR.detZ1   = 0;

end