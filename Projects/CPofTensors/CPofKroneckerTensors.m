%% Symbolic
B = sym('tB_%d_%d_%d', [2 2 2])
C = sym('tC_%d_%d_%d', [2 2 2])
A = superkron(B,C)

A = sym('tA_%d_%d_%d',[4 4 4])

%% Numerical
B = rand(2,2,2)
C = eyen(2,3)
A = superkron(C,B)

heig(B)
heig(C)
ea = heig(A)

A = superkron(B,C)
heig(A)

numSolns(4,3)

%% Random
B = rand(2,2,2);
C = rand(2,2,2);
A = superkron(B,C);

eb = heig(B)
ec = heig(C)
ea = heig(A)

ebc = kron(eb,ec)


unique([kron(eb, ec) kron(eb, kron(eb, ec)), kron(ec, kron(ec, eb))])

ea
eb
ec
ebc

%% Gleich example

N4 = eyen(4,3)
e4 = zeig(N4)
N2 = eyen(2,3)
e2 = zeig(N2)

e4

