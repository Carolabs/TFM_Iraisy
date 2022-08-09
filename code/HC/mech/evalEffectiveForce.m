function EffForce = evalEffectiveForce(MECH)

% Returns the effective force of the mechanical system along the direction
% of the hydraulic actuator

% In direct co-simulation, this function is used only to initialize the
% hydraulics

% ________________________________ Eval dynamics terms of mechanical system
MM          = evalMassMatrix(MECH);
fg          = evalGravForces(MECH);
fc          = evalVelDepForces(MECH);

% Eval Jacobian matrix of cylinder direction
[~, ~, J] = findSfromAngles(MECH);

meff        = (J*MM^(-1)*J')^(-1);
EffForce    = meff * J * MM^(-1) * (fg-fc);