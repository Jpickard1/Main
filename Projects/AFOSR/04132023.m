% 3-Uniform hyperstar on 5 nodes
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 13, 2023

clear
HG = getToyHG(5,3,'hyperstar');
[O, J] = HGObsvSym(HG)

simplify(O{3}(1:5,:))

rank(O{3}(:,1:2))

rank(O{3})

J{1}
J{2}
J{3}

rank(O{3}(1:4,:))
O{3}(1:4,:)

rank(O{2})
rank(O{2}(1:3,:))

r = [O{2}(1:3,:); O{2}(1:3,:)]
rank(r)
