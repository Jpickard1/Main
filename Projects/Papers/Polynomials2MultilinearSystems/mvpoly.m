classdef mvpoly
    %MVPOLY Multivariate Polynomial Systems
    %   This class represents and evaluates multivariate polynomial systems
    %
    % Auth: Joshua Pickard
    %       jpic@umich.edu
    % Date: August 5, 2023

    properties
        Am
        maxD
        nvars
    end

    methods
        function obj = mvpoly(NVA)
            %mvpoly Construct an instance of this class
            arguments
                NVA.Am
                NVA.maxD
                NVA.nvars
            end

            if ~isfield(NVA, 'Am');      NVA.Am      = []; end
            if ~isfield(NVA, 'maxD');    NVA.maxD    = []; end
            if ~isfield(NVA, 'nvars'); NVA.nvars     = []; end

            obj.Am      = NVA.Am;
            obj.maxD    = NVA.maxD;
            obj.nvars = NVA.nvars;
        end

        function Y = eval(obj, x)
            n = numel(x);
            nsv = numPolyTerms(n, obj.maxD);
            sv = zeros(nsv, 1);
            sv(nsv) = 1; nsv = nsv - 1;
            for i=1:obj.maxD
                sv(nsv-n^i+1:nsv) = KroneckerPower(x, i);
                nsv = nsv - n^i;
            end
            Y = obj.Am * sv;
        end

    end
end

