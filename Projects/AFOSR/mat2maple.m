function [outputArg1,outputArg2] = mat2maple(HG, filename)
%MAT2MAPLE This function generates a file for use by the MAPLE
%   implementation of Sedoglavic's algorithm.
%
%   INPUTS:
%       HG.........Hypergraph object
%       filename...name of output file
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 14, 2023
%%
filename = 't2.txt';
HG = HAT.uniformErdosRenyi(20,20,5);
disp(full(HG.IM));

% Extract uniform edge set;
ES = HAT.uniformEdgeSet(HG);
fileID = fopen(filename,'w');

output = {};
output{length(output) + 1} = 'infolevel[observabilityTest] := 1 :';
output{length(output) + 1} = 'infolevel[observabilitySymmetries] := 1 :';
output{length(output) + 1} = 'VectorField:= [';

% Write out equations
for e=1:size(ES,1)
    eqn = "    ";
    for vx=1:size(ES,2)
        if vx~=1
            eqn = eqn + "*";
        end
        eqn = eqn + "x" + string(ES(e,vx));
    end
    if e ~= size(ES,1)
        eqn = eqn + ",";
    end
    output{length(output) + 1} = eqn;
end
% End equations

output{length(output) + 1} = "]:"; 
output{length(output) + 1} = "OutputSystem := [";
output{length(output) + 1} = "x1";
output{length(output) + 1} = "] :";
output{length(output) + 1} = "OutputsVariables:= [y] :";
output{length(output) + 1} = "Inputs := [] :";
output{length(output) + 1} = "Parameters 	:= [] :";

% Write out variables
vars = "Variables 	:= [";
for i=1:size(HG.IM,2)
    if i~=1
        vars = vars + ",";
    end
    vars = vars + "x" + string(i);
end
vars = vars + "] :";
output{length(output) + 1} = vars;

output{length(output) + 1} = "libname := cat(currentdir()," + '"' + "/release" + '"' + "),libname :";
output{length(output) + 1} = "readlib(observabilityTest) :";
output{length(output) + 1} = "NonObservable := observabilityTest(	VectorField	,";
output{length(output) + 1} = "					Variables	,";
output{length(output) + 1} = "					OutputSystem	,";
output{length(output) + 1} = "					Parameters	,";
output{length(output) + 1} = "					Inputs			) :";
output{length(output) + 1} = "print(%) :";
output{length(output) + 1} = "GroupInfGen := observabilitySymmetries(	VectorField	,";
output{length(output) + 1} = "					Variables	,";
output{length(output) + 1} = "					OutputSystem	,";
output{length(output) + 1} = "					Parameters	,";
output{length(output) + 1} = "					Inputs		,";
output{length(output) + 1} = "					NonObservable		) :";
output{length(output) + 1} = "print(%) :";
output{length(output) + 1} = "quit :";

% Write to file
for i=1:length(output)
    fprintf(fileID, output{i} + "\n")
end

end

