
hm = 10e-3;     % Step mecánica
hh = 1e-4;      % Step hidraulica
st = 100;       % Step ratio
tfin = 10;      
coeffs = [];    % Coeficientes estimador
lim = 55;       % Variación máxima de la fuerza
ite = 10;       % Iteraciones

% Primera llamada al código de cosim sin aplicar corrección
[STORE_HYD, STORE_MECH] = main_Jacobi_fs(hm, hh, st, tfin, 0, coeffs, lim);

STORE_HYD_orig = STORE_HYD;
STORE_MECH_orig = STORE_MECH;

% ___________________________________________________ Monolithic
addpath('./monolithic')
mono_results = mainscript_monolithic(tfin, hm); % Esto da una solución de referencia
rmpath('./monolithic')

 
for i=1:ite

    % ____________________________________________________________ Estimation
    X = ones(1000, 3);
    X(:, 2) = STORE_HYD.F(1:1000)';
    X(:, 3) = STORE_HYD.F(2:1001)';
    Y = mono_results.F(1, 2:1001)';
    
    coeffs = estimator(X,Y); % Coeficientes del estimador usado para esta iteración
    
    % Repetir la llamada al código de cosim, pero corrigiendo con el nuevo estimador
    [STORE_HYD, STORE_MECH] = main_Jacobi_fs(hm, hh, st, tfin, 1, coeffs, lim);

end

% ____________________________________________________________ Plots
indexPlots = 0;

% Positions
indexPlots =  indexPlots + 1;
figure(indexPlots)
clf;
hold on
set(indexPlots, 'name', 'Actuator distance');
plot(mono_results.t, mono_results.pos(7,:));
plot(STORE_MECH.t, STORE_MECH.s, 'r--');
plot(STORE_MECH_orig.t, STORE_MECH_orig.s, 'g');
xlabel('Time (sec)');
ylabel('s (m)');
legend('mono', 'cosim', 'orig');
hold off

indexPlots =  indexPlots + 1;
figure(indexPlots)
clf;
hold on
set(indexPlots, 'name', 'Actuator velocity');
plot(mono_results.t, mono_results.vel(7,:));
plot(STORE_MECH.t, STORE_MECH.sd, 'r--');
plot(STORE_MECH_orig.t, STORE_MECH_orig.sd, 'g');
xlabel('Time (sec)');
ylabel('sd (m/s)');
legend('mono', 'cosim', 'orig');
hold off

% Hydraulic pressures
indexPlots =  indexPlots + 1;
figure(indexPlots)
clf;
hold on
set(indexPlots, 'name', 'Pressures');
plot(mono_results.t,mono_results.p(1,:)/1e6, 'r');
plot(STORE_HYD.t,STORE_HYD.p(1,:)/1e6, 'b');
plot(STORE_HYD_orig.t,STORE_HYD_orig.p(1,:)/1e6, 'g');
legend('mono','cosim', 'orig');
xlabel('Time (sec)');
ylabel('Pressure 1 (MPa)');
hold off

indexPlots =  indexPlots + 1;
figure(indexPlots)
clf;
hold on
set(indexPlots, 'name', 'Pressures');
plot(mono_results.t,mono_results.p(2,:)/1e6, 'r');
plot(STORE_HYD.t,STORE_HYD.p(2,:)/1e6, 'b');
plot(STORE_HYD_orig.t,STORE_HYD_orig.p(2,:)/1e6, 'g');
legend('mono','cosim', 'orig');
xlabel('Time (sec)');
ylabel('Pressure 2 (MPa)');
hold off

% Coupling forces
indexPlots =  indexPlots + 1;
figure(indexPlots)
clf;
hold on
set(indexPlots, 'name', 'Actuator force');
plot(mono_results.t, mono_results.F);
plot(STORE_HYD.t, STORE_HYD.F, 'r--');
plot(STORE_HYD_orig.t, STORE_HYD_orig.F, 'g');
xlabel('Time (sec)');
ylabel('Force (N)');
legend('hyd', 'mech', 'orig');
hold off