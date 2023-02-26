function [outputArg1,outputArg2] = mat2maple(HG, outputVars, filename)
%MAT2MAPLE This function generates a file for use by the MAPLE
%   implementation of Sedoglavic's algorithm.
%
%   INPUTS:
%       HG.........Hypergraph object
%       filename...name of output file
%
%   EXAMPLE
%{
    V = 5;k=3;
    HG = hyperring(V,k);
    outputVars = [1 2];
    filename = "t7.mpl"
    mat2maple(HG, outputVars, filename)
%}
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 14, 2023
%%
% filename = 't2.mpl';
% HG = HAT.uniformErdosRenyi(20,20,5);
disp(full(HG.IM));

% Extract uniform edge set;
% ES = HAT.uniformEdgeSet(HG);
output = {};
output{length(output) + 1} = 'infolevel[observabilityTest] := 1 :';
output{length(output) + 1} = 'infolevel[observabilitySymmetries] := 1 :';
output{length(output) + 1} = 'VectorField:= [';

% Write out equations
for vx=1:size(HG.IM,1)
    e = find(HG.IM(vx,:) == 1);
    eqn = "    ";
    for ei=1:length(e)
        if ~strcmp(eqn, "    ")
            eqn = eqn + "+";
        end
        vxc = find(HG.IM(:,e(ei)) == 1);
        vxc = vxc(vxc ~= vx);
        for i=1:length(vxc)
            if i ~= 1
                eqn = eqn + "*";
            end
            eqn = eqn + "x" + string(vxc(i));
        end
    end
    if vx ~= size(HG.IM,1)
        eqn = eqn + ",";
    end
    output{length(output) + 1} = eqn;
end
% End equations

output{length(output) + 1} = "]:"; 
output{length(output) + 1} = "OutputSystem := [";
for i=1:length(outputVars)
    if i ~= length(outputVars)
        output{length(output) + 1} = "    x" + string(outputVars(i)) + ",";
    else
        output{length(output) + 1} = "    x" + string(outputVars(i));
    end
end
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
output{length(output) + 1} = "print(#) :";
output{length(output) + 1} = "GroupInfGen := observabilitySymmetries(	VectorField	,";
output{length(output) + 1} = "					Variables	,";
output{length(output) + 1} = "					OutputSystem	,";
output{length(output) + 1} = "					Parameters	,";
output{length(output) + 1} = "					Inputs		,";
output{length(output) + 1} = "					NonObservable		) :";
output{length(output) + 1} = "print(%) :";
output{length(output) + 1} = "quit :";

% Write to file
fileID = fopen(filename,'w');
for i=1:length(output)
    fprintf(fileID, output{i} + "\n");
end

end

%{
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
%}
