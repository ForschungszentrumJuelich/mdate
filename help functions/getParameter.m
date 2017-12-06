function [ varargout ] = getParameter(model, p, lambda)
%% to do
% create following calling structure:
    % function [ varargout ] = getParameter(varargin)
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% created 10.08.2011 by Jonas Bühler, Nico Becker, Eric von Lieres
%
% This function returns a structure containing all parameters for the
% respective model modelNo. There are two ways to call getParameter:
%
% calling interface:
%
% %   P = getParameter(model, p, lambda);
% %       integrate model parameters in p into the parameter structure P and
% %       return P
% %
% %   [pName pStart LB UB UB_MS] = getParameter(model, 'pInfo')
% %       THIS FUNCTION IS MORE A LOOK UP TABLE THAN A FUNCTION!
% %       returns parameter names of model modelNo in cell array P. LB and UB
% %       contain upper and lower boundaries for fitting; these values might
% %       not be unversally valid but can give a hint of the orders of
% %       magnitude of parameter ranges. pStart might be a helpful starting
% %       parameter for fitting but is not unversally valid also.
% %       
% %       pName: Name of each model parameter, e.g. u_1 (velocity in compartment 1)
% %       pStart: suggested starting value for optimization
% %       LB: lower boundaries for optimization
% %       UB: upper boundaries for optimization
% %       UB_MS: Upper boundaries for starting values in multistart optimization
%

%  annotation from 20120704
%     - extend getParameter such to allow the calculation of a model given
%       by the connectivity matrix
%     - alternatively create a function that creates a typical model
%       structure out of a given connectivity matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setting model meta information N, BC, IC, ...
while 1
    
    % split model description entry
    [~,~,~,~,~,~,splitstring] = regexp(model{1}, '_');

    % get name of model property/parameter
    pname = splitstring{1};
    
    % set values
    switch pname
        case 'N'
            N = str2double(splitstring{2});
            sscanf(model{1}, 'N_%d'); %% loeschen??? tut nichts
            
            %% model definition
            A = zeros(N,N);
            u = zeros(1,N);
            d = zeros(1,N);
            %%
            
            %% default initial and boundary function pointer (zero)
            inletFun   = cell(N,1);
            initialFun = cell(N,1);
            for i = 1:N
%                 inletFun{i}   = @inlet_zero;
%                 initialFun{i} = @initial_zero;
                inletFun{i}   = 'zero';
                initialFun{i} = 'zero';
            end
            %%
            
            %% Spatial Gaussian initial condition
            x0     = zeros(N,1);
            xSigma = zeros(N,1);
            %%
            
            %% Temporal Gaussian boundary condition
            t0     = zeros(N,1);
            tSigma = zeros(N,1);
            %%
            
            %% Temporal square boundary condition
            t1     = zeros(N,1);
            t2     = zeros(N,1);
            %%
            
            %% Spatial square boundary condition
            x1     = zeros(N,1);
            x2     = zeros(N,1);
            %%
            
            %% Scaling of initial and boundary condition
            IC0 = zeros(N,1);
            BC0 = zeros(N,1);
            %%
            
            %% steepness factor of smoothed edge
            edgeSteepNess = zeros(N,1);
            %%
            
        case 'BC'
            i = str2double(splitstring{2});
            fun = splitstring{3};
            BC0(i) = 1;
            switch fun
                case 'data'
                    inletFun{i}   = 'data';
                case 'gauss'
                    inletFun{i}   = 'gauss';
                case 'smoothEdge'
                    inletFun{i}   = 'smoothEdge';
                case 'square'
                    inletFun{i}   = 'square';
                case 'zero'
                    inletFun{i}   = 'zero';
            end
        case 'IC'
            i = str2double(splitstring{2});
            fun = splitstring{3};
            IC0(i) = 1;
            switch fun
                case 'data'
                    initialFun{i}   = 'data';
                case 'gauss'
                    initialFun{i}   = 'gauss';
                case 'square'
                    initialFun{i}   = 'square';
                case 'smoothEdge'
                    initialFun{i}   = 'smoothEdge';
                case 'zero'
                    initialFun{i}   = 'zero';
            end
        otherwise
            % all meta information about model have to be the first
            % entries. at this point the model description should be of the
            % same size as p and only contain the respective parameter
            % names
            if (length(model) ~= length(p))  &&   ~ischar(p)   && ~isempty(p)
                error('parameter list in model and length(p) are unequal!');
            end
            break
    end
    
    % remove first entry from model description after reading it
    model(1) = [];
    
end
%%

%% define number of model parameters
np = length(model);
%%

%% function call #2 "pInfo", setting LB UB
% return parameter names, suggested starting values and upper and lower boundaries
if strcmp(p, 'pInfo')
    %% parameter names
    pStart = zeros(np, 1);
    LB     = zeros(np, 1);
    UB     = zeros(np, 1);
    LB_MS  = zeros(np, 1);
    UB_MS  = zeros(np, 1);
    
    %% hardcoded starting and boundary values 
    start_u      = 1.0;       LB_u      = 1e-1;     LB_MS_u      = 0.5;     UB_u      = 10;    UB_MS_u        = 5;
    start_d      = 1e-10;     LB_d      = 1e-20;    LB_MS_d      = 1e-10;   UB_d      = 1;      UB_MS_d        = 1;
    start_a      = 0.1;       LB_a      = 1e-8;     LB_MS_a      = LB_a;    UB_a      = 2;    UB_MS_a        = 1;
    start_BC0    = 1;         LB_BC0    = 1e-8;     LB_MS_BC0    = LB_BC0;  UB_BC0    = 1e3;    UB_MS_BC0      = 1;
    start_IC0    = 1;         LB_IC0    = 1e-8;     LB_MS_IC0    = 0.5;     UB_IC0    = 1e3;    UB_MS_IC0      = 1;
    start_x0     = 40;        LB_x0     = 10;       LB_MS_x0     = LB_x0;   UB_x0     = 300;    UB_MS_x0       = 200;
    start_xSigma = 5;         LB_xSigma = 3;        LB_MS_xSigma = 0.1;     UB_xSigma = 100;    UB_MS_xSigma   = 30;
    start_t0     = 50;        LB_t0     = 1;        LB_MS_t0     = 1;       UB_t0     = 1e3;    UB_MS_t0       = 100;
    start_tSigma = 15;        LB_tSigma = 1;        LB_MS_tSigma = 1;       UB_tSigma = 200;    UB_MS_tSigma   = 30;
    start_t1     = 20;        LB_t1     = 0;        LB_MS_t1     = 0.5;     UB_t1     = 1e3;    UB_MS_t1       = 0;
    start_t2     = 30;        LB_t2     = 0;        LB_MS_t2     = 0.5;     UB_t2     = 1e3;    UB_MS_t2       = 0;
    start_x1     = 50;        LB_x1     = 0;        LB_MS_x1     = 0.5;     UB_x1     = 1e3;    UB_MS_x1       = 0;
    start_x2     = 60;        LB_x2     = 0;        LB_MS_x2     = 0.5;     UB_x2     = 1e3;    UB_MS_x2       = 0 ;
    start_eSN    = 1;         LB_eSN    = 0.1;      LB_MS_eSN    = 0.1;     UB_eSN    = 10;    UB_MS_eSN       = 10 ;
    start_ampl   = 1;         LB_ampl   = .5;       LB_MS_ampl   = 0.95;    UB_ampl   = 2.0;    UB_MS_ampl     = 1.05;
    %%
    
    % set output parameter
    for k = 1:np;
        % get parameter name
        pname = strtok(model{k}, '_');
        switch pname
            case 'u'
                pStart(k) = start_u;
                LB(k)     =    LB_u;
                UB(k)     =    UB_u;
                LB_MS(k)  = LB_MS_u;
                UB_MS(k)  = UB_MS_u;
                
            case 'd'
                pStart(k) = start_d;
                LB(k)     =    LB_d;
                UB(k)     =    UB_d;
                LB_MS(k)  = LB_MS_d;
                UB_MS(k)  = UB_MS_d;
            
            case 'a'
                pStart(k) = start_a;
                LB(k)     =    LB_a;
                UB(k)     =    UB_a;
                LB_MS(k)  = LB_MS_a;
                UB_MS(k)  = UB_MS_a;
                
            case 'BC0'
                pStart(k) = start_BC0;
                LB(k)     =    LB_BC0;
                UB(k)     =    UB_BC0;
                LB_MS(k)  = LB_MS_BC0;
                UB_MS(k)  = UB_MS_BC0;

            case 'IC0'
                pStart(k) = start_IC0;
                LB(k)     =    LB_IC0;
                UB(k)     =    UB_IC0;
                LB_MS(k)  = LB_MS_IC0;
                UB_MS(k)  = UB_MS_IC0;
            
            case 'x0'
                pStart(k) = start_x0;
                LB(k)     =    LB_x0;
                UB(k)     =    UB_x0;
                LB_MS(k)  = LB_MS_x0;
                UB_MS(k)  = UB_MS_x0;

            case 'xSigma'
                pStart(k) = start_xSigma;
                LB(k)     =    LB_xSigma;
                UB(k)     =    UB_xSigma;
                LB_MS(k)  = LB_MS_xSigma;
                UB_MS(k)  = UB_MS_xSigma;

            case 't0'
                pStart(k) = start_t0;
                LB(k)     =    LB_t0;
                UB(k)     =    UB_t0;
                LB_MS(k)  = LB_MS_t0;
                UB_MS(k)  = UB_MS_t0;

            case 'tSigma'
                pStart(k) = start_tSigma;
                LB(k)     =    LB_tSigma;
                UB(k)     =    UB_tSigma;
                LB_MS(k)  = LB_MS_tSigma;
                UB_MS(k)  = UB_MS_tSigma;
                
            case 't1'
                pStart(k) = start_t1;
                LB(k)     =    LB_t1;
                UB(k)     =    UB_t1;
                LB_MS(k)  = LB_MS_t1;
                UB_MS(k)  = UB_MS_t1;
                
            case 't2'
                pStart(k) = start_t2;
                LB(k)     =    LB_t2;
                UB(k)     =    UB_t2;
                LB_MS(k)  = LB_MS_t2;
                UB_MS(k)  = UB_MS_t2;
                
            case 'x1'
                pStart(k) = start_x1;
                LB(k)     =    LB_x1;
                UB(k)     =    UB_x1;
                LB_MS(k)  = LB_MS_x1;
                UB_MS(k)  = UB_MS_x1;
                
            case 'x2'
                pStart(k) = start_x2;
                LB(k)     =    LB_x2;
                UB(k)     =    UB_x2;
                LB_MS(k)  = LB_MS_x2;
                UB_MS(k)  = UB_MS_x2;
                
            case 'edgeSteepNess'
                pStart(k) = start_eSN;
                LB(k)     =    LB_eSN;
                UB(k)     =    UB_eSN;
                LB_MS(k)  = LB_MS_eSN;
                UB_MS(k)  = UB_MS_eSN;
                
            case 'ampl'
                pStart(k) = start_ampl;
                LB(k)     =    LB_ampl;
                UB(k)     =    UB_ampl;
                LB_MS(k)  = LB_MS_ampl;
                UB_MS(k)  = UB_MS_ampl;
                
            case 'fourier'
                pStart = [];
                LB = [LB_u LB_a LB_a LB_a LB_x0 LB_xSigma LB_ampl]';
                UB = [UB_u UB_a UB_a UB_a UB_x0 UB_xSigma UB_ampl]';
                LB_MS = [LB_MS_u LB_MS_a LB_MS_a LB_MS_a LB_MS_x0 LB_MS_xSigma LB_MS_ampl]';
                UB_MS = [UB_MS_u UB_MS_a UB_MS_a UB_MS_a UB_MS_x0 UB_MS_xSigma UB_MS_ampl]';
                
                
            otherwise
                error( ['Unknown parameter: ' pname] )
        end
    end
    
    if nargout == 5
        error('Using probably old interface for [model, pStart, LB, UB, LB_MS, UB_MS] = getParameter(model, ''pInfo'')... There should be 6 lhs parameters.');
    end
    
    varargout{1} = model;
    varargout{2} = pStart;
    varargout{3} = LB;
    varargout{4} = UB;
    varargout{5} = LB_MS;
    varargout{6} = UB_MS;
    %%
    return
end
%% eo "pInfo" section

%% set default for amplitude - efficient solver paper only
amplitude = 1;
%%

%% function call #1  P = getParameter(model, p, lambda);
% setting model parameter values
for k = 1:np
    
    % get parameter name
    pname = strtok(model{k}, '_');
    parameter = p(k);
    switch pname
        case 'u'
            i = sscanf(model{k}, 'u_%d');
            u(i) = parameter;
        case 'd'
            i = sscanf(model{k}, 'd_%d');
            d(i) = parameter;
        case 'a'
            idx = sscanf(model{k}, 'a_%d_%d');
            i = idx(1);
            j = idx(2);
            A(i,i) = A(i,i) - parameter; % subtract from diagonal
            A(j,i) = A(j,i) + parameter; % add to respective compartment
        case 'BC0'
            i = sscanf(model{k}, 'BC0_%d');
            BC0(i) = parameter;
        case 'IC0'
            i = sscanf(model{k}, 'IC0_%d');
            IC0(i) = parameter;
        case 'x0'
            i = sscanf(model{k}, 'x0_%d');
            x0(i) = parameter;
        case 'xSigma'
            i = sscanf(model{k}, 'xSigma_%d');
            xSigma(i) = parameter;
        case 't0'
            i = sscanf(model{k}, 't0_%d');
            t0(i) = parameter;
        case 'tSigma'
            i = sscanf(model{k}, 'tSigma_%d');
            tSigma(i) = parameter;
        case 't1' % begin of square function
            i = sscanf(model{k}, 't1_%d');
            t1(i) = parameter;
        case 't2' % end of square function
            i = sscanf(model{k}, 't2_%d');
            t2(i) = parameter;
        case 'x1' % begin of square function
            i = sscanf(model{k}, 'x1_%d');
            x1(i) = parameter;
        case 'x2' % end of square function
            i = sscanf(model{k}, 'x2_%d');
            x2(i) = parameter;
        case 'edgeSteepNess'
            i = sscanf(model{k}, 'edgeSteepNess_%d');
            edgeSteepNess(i) = parameter;
        case 'ampl'
            i = sscanf(model{k}, 'ampl_%d');
            if isempty(i)
                amplitude = parameter;
            else
                amplitude(i) = parameter;
            end
        otherwise
            error(['unknown parameter name: ' pname]);
    end
end
%%

%% consistency checks for boundary and initital conditions
for i=2:N
    % if no parameter for BC/IC in compartment i exist, use BC/IC parameter
    % from previous compartment NEEDS TO BE EXTENDEND FOR SQUARE FUNCTION
    % in Intital Condition
    if BC0(i) > 0 && tSigma(i) == 0
        tSigma(i) = tSigma(i-1);
        t0(i)     = t0(i-1);
    end
    if BC0(i) > 0 && t1(i) == 0
        t1(i) = t1(i-1);
        t2(i) = t2(i-1);
    end
    if IC0(i) > 0 && xSigma(i) == 0
        xSigma(i) = xSigma(i-1);
        x0(i)     = x0(i-1);
    end
    if IC0(i) > 0 && x1(i) == 0
        x1(i) = x1(i-1);
        x2(i) = x2(i-1);
    end
end
%%

%% subtract decay rate
A = A - lambda*eye(N);
%%

P.u          = u;
P.d          = d;
P.N          = N;
P.A          = A;
P.inletFun   = inletFun;
P.initialFun = initialFun;
P.IC0        = IC0;
P.BC0        = BC0;
P.x0         = x0;
P.xSigma     = xSigma;
P.t0         = t0;
P.tSigma     = tSigma;
P.t1         = t1;
P.t2         = t2;
P.x1         = x1;
P.x2         = x2;
P.edgeSteepNess = edgeSteepNess;
P.amplitude  = amplitude;
P.lambda     = lambda;



varargout{1} = P;

