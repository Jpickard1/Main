classdef multilinearSystem
    %MULTILINEARSYSTEM MULTILINEAR TIME INVARIANT SYSTEM
    %   This class represents and evaluates multilinear systems
    %
    % Auth: Joshua Pickard
    %       jpic@umich.edu
    % Date: August 5, 2023

    properties
        A
        nvars
    end
    
    methods
        function obj = multilinearSystem(NVA)
            %MULTILINEARSYSTEM Construct an instance of this class
            arguments
                NVA.A
                NVA.poly
            end

            if isfield(NVA, 'poly')
                obj.A = poly2tensor(NVA.poly);
            elseif isfield(NVA, 'A')
                obj.A = NVA.A;
            else
                error('Either a polynomial or a tensor is required for construction.');
            end
        end
        
        function Y = eval(obj, x)
            x = [x; 1];
            Y = obj.A;
            while numel(double(Y)) > numel(x)
                Y = ttv(Y, x, ndims(Y));
            end
            Y = double(Y);
            Y = Y(1:numel(Y)-1);
        end
    end
end

