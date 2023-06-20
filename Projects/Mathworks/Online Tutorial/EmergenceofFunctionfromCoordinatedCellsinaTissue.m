% EQ1A = sbioabstractkineticlaw('EQ1A', 'a * (x - x0)')
% set(EQ1A, 'SpeciesVariables', {'x'});
% set(EQ1A, 'ParameterVariables', {'a', 'x0'});
% sbioaddtolibrary(EQ1A);
% EQ1B = sbioabstractkineticlaw('EQ1B', 'beta * (x - y)')
% set(EQ1B, 'SpeciesVariables', {'x', 'y'});
% set(EQ1B, 'ParameterVariables', {'beta'});
% sbioaddtolibrary(EQ1B);

clear all; close all; clc

%% Equations 2 from the paper
m1 = sbiomodel('simpleExample');                  % Create a SimBiology Model
addspecies(m1, 'x', 'InitialAmount', 10)
addspecies(m1, 'y', 'InitialAmount', 1)

r1 = addreaction(m1, 'x -> x');
r2 = addreaction(m1, 'y -> y');
r3 = addreaction(m1, 'x -> y');
r4 = addreaction(m1, 'y -> x');

klaw1 = addkineticlaw(r1, 'EQ1A', 'ParameterVariableNames', {'a','x0'}, 'SpeciesVariableNames', {'x'})
klaw2 = addkineticlaw(r2, 'EQ1A', 'ParameterVariableNames', {'a','y0'}, 'SpeciesVariableNames', {'y'})
klaw3 = addkineticlaw(r3, 'EQ1B', 'ParameterVariableNames', {'beta'}, 'SpeciesVariableNames', {'x','y'})
klaw4 = addkineticlaw(r4, 'EQ1B', 'ParameterVariableNames', {'beta'}, 'SpeciesVariableNames', {'y','x'})

a = 0.5;
b = 0.5;
x0 = 2;
y0 = 2;
beta = 1;

addparameter(klaw1,'a',a);
addparameter(klaw1,'x0',x0);
addparameter(klaw2,'a',a);
addparameter(klaw2,'y0',y0);
addparameter(klaw3,'beta',beta);
addparameter(klaw4,'beta',beta);

% Simulate model:
sd = sbiosimulate(m1);
% Plot results:
sbioplot(sd);
% Open the plot in model viewer
simbiology(m1)

%% Copy compartments
clear all; close all; clc
m2 = sbiomodel('duplicateCompartments');                  % Create a SimBiology Model
c1 = addcompartment(m2, 'c1');
% c2 = addcompartment(m2, 'c2')

addspecies(c1, 'x', 'InitialAmount', 10);
addspecies(c1, 'y', 'InitialAmount', 1);

r1 = addreaction(m2, 'c1.x -> c1.x');
r2 = addreaction(m2, 'c1.y -> c1.y');
r3 = addreaction(m2, 'c1.x -> c1.y');
r4 = addreaction(m2, 'c1.y -> c1.x');

klaw1 = addkineticlaw(r1, 'EQ1A', 'ParameterVariableNames', {'a','x0'}, 'SpeciesVariableNames', {'x'});
klaw2 = addkineticlaw(r2, 'EQ1A', 'ParameterVariableNames', {'a','y0'}, 'SpeciesVariableNames', {'y'});
klaw3 = addkineticlaw(r3, 'EQ1B', 'ParameterVariableNames', {'beta'}, 'SpeciesVariableNames', {'x','y'});
klaw4 = addkineticlaw(r4, 'EQ1B', 'ParameterVariableNames', {'beta'}, 'SpeciesVariableNames', {'y','x'});

a = 0.5;
b = 0.5;
x0 = 2;
y0 = 2;
beta = 1;

addparameter(klaw1,'a',a);
addparameter(klaw1,'x0',x0);
addparameter(klaw2,'a',a);
addparameter(klaw2,'y0',y0);
addparameter(klaw3,'beta',beta);
addparameter(klaw4,'beta',beta);

r11 = copyobj(r1, m2)
r21 = copyobj(r2, m2)
r31 = copyobj(r3, m2)
r41 = copyobj(r4, m2)

c2 = copyobj(c1, m2);
m2.addcompartment(c2)

simbiology(m2)


c1 = m1.Compartments(1)
m1.Compartments(2) = c1;
simbiology(m1)


%% 
m2.Species
% Simulate model:
sd = sbiosimulate(m2);
% Plot results:
sbioplot(sd);

m1.Compartments