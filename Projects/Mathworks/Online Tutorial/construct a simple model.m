% https://www.mathworks.com/help/simbio/gs/construct-a-simple-model.html

m1 = sbiomodel('simpleModel');                  % Create a SimBiology Model
r1 = addreaction(m1,'A -> B');                  % Add a reaction to  
m1.species                                      % View species in m1
m1.species(1).InitialAmount = 10;               % Set initial conditions
kineticLaw = addkineticlaw(r1,'MassAction');    % Define form of r1
p1 = addparameter(kineticLaw,'k',0.5);          % Define parameter of r1
kineticLaw.ParameterVariableNames = 'k';        % Set param name of r1
sd = sbiosimulate(m1);                          % Simulate the model
sbioplot(sd);
