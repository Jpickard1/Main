% https://drive.google.com/file/d/1zH4n_rsLeF3x-HVCrkujALil7272Tg3P/view
% https://drive.google.com/file/d/1VS-GnZR0djkldLk9yJ4V3OfpynBz1qVl/view

m1 = sbiomodel('simpleExample');                  % Create a SimBiology Model
addspecies(m1, 'x1', 'InitialAmount', 10)
addspecies(m1, 'y1', 'InitialAmount', 1)
addspecies(m1, 'x2', 'InitialAmount', 10)
addspecies(m1, 'y2', 'InitialAmount', 1)


r1 = addreaction(m1, 'x -> x');
r2 = addreaction(m1, 'y -> y');
r3 = addreaction(m1, 'x -> y');
r4 = addreaction(m1, 'y -> x');

matrixToEquations()

help equationsToMatrix

