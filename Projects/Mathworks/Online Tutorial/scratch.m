
%% 
p1 = addparameter(klaw1,'a',a);                 % Define parameter of r1
p2 = addparameter(klaw1,'x0',x0);               % Define parameter of r1
klaw1.ParameterVariableNames = {'a', 'x0'};     % Set param name of r1
klaw1.SpeciesVariables = {'x'};

p3 = addparameter(klaw2,'a',b);                 % Define parameter of r1
p4 = addparameter(klaw2,'x0',y0);               % Define parameter of r1
klaw2.ParameterVariableNames = {'a', 'x0'};     % Set param name of r1
klaw2.SpeciesVariables = {'y'};

p5 = addparameter(klaw3,'beta',beta);              % Define parameter of r1
klaw3.ParameterVariableNames = {'beta'};  % Set param name of r1
klaw3.SpeciesVariables = {'x','y'};

p6 = addparameter(klaw4,'beta',beta);              % Define parameter of r1
klaw4.ParameterVariableNames = {'beta'};  % Set param name of r1
klaw4.SpeciesVariables = {'x','y'};


m1.species(1).InitialAmount = 10;               % Set initial conditions
m1.species(2).InitialAmount = 5;               % Set initial conditions
sd = sbiosimulate(m1);                          % Simulate the model
sbioplot(sd);

simbiology(m1);

%% Online help
% Create model and add species:
M = sbiomodel ('CI model');
S = addspecies(M, 'Substrate', 'InitialAmount', 10);
E = addspecies(M, 'Enzyme', 'InitialAmount', 1);
P = addspecies(M, 'Product', 'InitialAmount', 0);
I = addspecies(M, 'Inhibitor', 'InitialAmount', 10);
% Add reaction to model:
R = addreaction(M, 'Substrate + Enzyme <-> Enzyme + Product');
% Add kinetic law:
KineticLaw = addkineticlaw(R, 'Competitive-Inhibition', ...
    'ParameterVariableNames', {'Km', 'Ki', 'Vm'}, ...
    'SpeciesVariableNames', {'Substrate', 'Inhibitor'});  
% Add parameters, whose names have to match the ones specified in
%  ParameterVariableNames, but they do not necessarily have to be
%  Km, Ki, or Vm in general.
Km = addparameter(KineticLaw, 'Km', 20);
Ki = addparameter(KineticLaw, 'Ki', 6);
Vm = addparameter(KineticLaw, 'Vm', 2);
% Simulate model:
sd = sbiosimulate(M);
% Plot results:
sbioplot(sd);


%% Intial try
clear
m1 = sbiomodel('simpleExample');                  % Create a SimBiology Model
r1 = addreaction(m1, 'x -> x');
r2 = addreaction(m1, 'y -> y');
r3 = addreaction(m1, 'x -> y');
r4 = addreaction(m1, 'y -> x');
m1.Species

EQ1A = sbioabstractkineticlaw('EQ1A', 'a * (x - x0)')
set(EQ1A, 'SpeciesVariables', {'x'});
set(EQ1A, 'ParameterVariables', {'a', 'x0'});
sbioaddtolibrary(EQ1A);

klaw1 = addkineticlaw(r1, 'EQ1A')
klaw2 = addkineticlaw(r2, 'EQ1A')
klaw3 = addkineticlaw(r3, 'EQ1A')
klaw4 = addkineticlaw(r4, 'EQ1A')

a = 0.5;
b = 0.5;
x0 = 2;
y0 = 2;
beta = 1;

p1 = addparameter(klaw1,'a',a);          % Define parameter of r1
p2 = addparameter(klaw1,'x0',x0);          % Define parameter of r1
klaw1.ParameterVariableNames = {'a', 'x0'};        % Set param name of r1
% klaw1.ParameterVariableNames = 'x0';        % Set param name of r1

p3 = addparameter(klaw2,'a',a);          % Define parameter of r1
p4 = addparameter(klaw2,'x0',y0);          % Define parameter of r1
klaw2.ParameterVariableNames = {'b', 'y0'};        % Set param name of r1
% klaw2.ParameterVariableNames = 'y0';        % Set param name of r1

p5 = addparameter(klaw3,'a',beta);          % Define parameter of r1
p6 = addparameter(klaw3,'x0',x0);          % Define parameter of r1
klaw3.ParameterVariableNames = {'beta', 'x0'};        % Set param name of r1
% klaw3.ParameterVariableNames = 'x0';        % Set param name of r1

p7 = addparameter(klaw4,'a',beta);          % Define parameter of r1
p8 = addparameter(klaw4,'x0',y0);          % Define parameter of r1
klaw3.ParameterVariableNames = {'beta', 'y0'};        % Set param name of r1
% klaw3.ParameterVariableNames = 'y0';        % Set param name of r1

simbiology(m1)
sbioremovefromlibrary('kineticlaw','EQ1A')

clear