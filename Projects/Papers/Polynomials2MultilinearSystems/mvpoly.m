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
        type
        args
    end

    methods
        function obj = mvpoly(NVA)
            %mvpoly Construct an instance of this class
            arguments
                NVA.Am      % matrix structure
                NVA.maxD    % maximum degree
                NVA.nvars   % number of variables
                NVA.type    % types of systems: lorenz, van der pol
                NVA.sigma   % lorenz
                NVA.rho     % lorenz
                NVA.beta    % lorenz, SIS
                NVA.mu      % van der pol
                NVA.gamma   % SIS
                NVA.stoch   % SIS
            end
            if isfield(NVA, 'type')
                if strcmp(NVA.type, 'lorenz')
                    NVA.Am    = lorenz(NVA.sigma, NVA.rho, NVA.beta);
                    NVA.maxD  = 2;
                    NVA.nvars = 3;
                    NVA.args  = {NVA.sigma, NVA.rho, NVA.beta};
                elseif strcmp(NVA.type, 'van der pol')
                    NVA.Am    = vanderpol(NVA.mu);
                    NVA.maxD  = 3;
                    NVA.nvars = 2;
                    NVA.args  = {NVA.mu};
                elseif strcmp(NVA.type, 'SIS')
                    NVA.Am    = SIS(NVA.beta, NVA.gamma, NVA.stoch);
                    NVA.maxD  = 2;
                    NVA.nvars = 2;
                    NVA.args  = {NVA.beta, NVA.gamma, NVA.stoch};
                elseif strcmp(NVA.type, 'JC2')
                    NVA.Am    = JCex2();
                    NVA.maxD  = 3;
                    NVA.nvars = 2;
                    NVA.args  = {};
                elseif strcmp(NVA.type, 'JB22')
                    NVA.Am    = JBex22();
                    NVA.maxD  = 2;
                    NVA.nvars = 2;
                    NVA.args  = {};
                else
                    error('invalid type');
                end
            end
            if ~isfield(NVA, 'Am');    NVA.Am      = [];    end
            if ~isfield(NVA, 'maxD');  NVA.maxD    = [];    end
            if ~isfield(NVA, 'nvars'); NVA.nvars   = [];    end
            if ~isfield(NVA, 'type');  NVA.type    = 'N/A'; end
            if ~isfield(NVA, 'args');  NVA.args    = {};    end

            obj.Am      = NVA.Am;
            obj.maxD    = NVA.maxD;
            obj.nvars   = NVA.nvars;
            obj.type    = NVA.type;
            obj.args    = NVA.args;
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

        function str = title(obj)
            if strcmp(obj.type, 'lorenz')
                str = "Lorenz ($\sigma=" + string(obj.args{1}) + ",\rho=" + string(obj.args{2}) + ",\beta=" + string(obj.args{3}) + "$)";
            elseif strcmp(obj.type, 'van der pol')
                str = "Van Der Pol ($\mu=" + string(obj.args{1}) + "$)";
            else
                str = "";
            end
        end

    end
end

