function MGR = evalManager(MGR, y1, y2, t)

% y1: SS1 outputs: s, sdot
% y2: SS2 outputs: force

% _________________________________________________________________________
%                                                      Store system outputs
MGR.storeIdx = MGR.storeIdx + 1;

MGR.STORE.t(MGR.storeIdx, 1)       = t;
MGR.STORE.s(MGR.storeIdx, 1)       = y1(1);
MGR.STORE.sd(MGR.storeIdx, 1)      = y1(2);
MGR.STORE.f(MGR.storeIdx, 1)       = y2(1);
MGR.STORE.fCorr(MGR.storeIdx, 1)   = MGR.fCorr; 
MGR.STORE.dP(MGR.storeIdx, 1)      = MGR.deltaP;
MGR.STORE.dE(MGR.storeIdx, 1)      = MGR.deltaE;
MGR.STORE.Ek(MGR.storeIdx, 1)      = MGR.Ek;
MGR.STORE.Ebalance(MGR.storeIdx, 1) = MGR.Ebalance;


% _________________________________________________________________________
%                                                       Evaluate indicators

% Actuator rate
xd  = y1(2);

% Using ZOH
deltaP1     = MGR.uTilde1 * y1(2);  % Force times velocity
deltaP2     = MGR.uTilde2 * y2(1);  % Velocity times force
MGR.deltaP  = -(deltaP1 + deltaP2);
MGR.deltaE  = MGR.deltaP * MGR.H;
MGR.Ek      = MGR.Ek + MGR.deltaE;

% Actual work of the corrective force
if (MGR.PR.correctEnergy == 1)
    MGR.dWfCorr = MGR.fCorr*xd * MGR.H;
else
    MGR.dWfCorr = 0;
end

% Energy error/2 in current step plus work of correction force
MGR.Ebalance = MGR.Ebalance + MGR.deltaE/2.0 + MGR.dWfCorr;

% _________________________________________________________________________
%                                                            Correct energy
if (MGR.PR.correctEnergy == 1)
    
    % Eliminate dEk/2
    
    % Required correction force
    if ((abs(MGR.deltaE) > MGR.PR.dEThold) && (abs(xd)>=MGR.PR.xdThold))
        fCorr 	= (MGR.PR.kp * MGR.deltaE + MGR.PR.kint*MGR.Ebalance) / (xd * MGR.H);
    else
        fCorr   = 0.0;
    end
    
    % Cap the force (max. correction force = corrFactor * coupling force)
    fCorrCap    = MGR.PR.fCapRel * abs(y2(1));
    
    if (fCorr >= fCorrCap) 
        %disp('Capping')
        fCorr = fCorrCap; 
    elseif (fCorr <= -fCorrCap) 
        %disp('Capping')
        fCorr = -fCorrCap;
    end
    
    % Change sign to comply with signs of coupling force
    MGR.fCorr   = -fCorr;
        
end

% _________________________________________________________________________
%                                               Update inputs to subsystems

% Outputs actually exchanged between subsystems
% These do not include the corrections
yEx     = [y1(1);y1(2);y2(1)];
MGR.uEx = MGR.L * yEx; 

% Update extrapolated inputs
% These do not include corrections and are used to evaluate the indicator
MGR.uTilde1 = MGR.uEx(1);
MGR.uTilde2 = MGR.uEx(3);

% Apply corrections if necessary
if (MGR.PR.correctEnergy == 1)
    yEx     = [y1(1);y1(2); y2(1) + MGR.fCorr];
    
    % Applied force cannot exceed absolute limit
    if ( abs(y2(1) + MGR.fCorr) > MGR.PR.fCapAbs)
        yEx(3) = (y2(1) + MGR.fCorr) * MGR.PR.fCapAbs / abs(y2(1) + MGR.fCorr);
    end
    
    MGR.uEx = MGR.L * yEx; 
end



end