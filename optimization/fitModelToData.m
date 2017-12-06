function [pFit, chiSqr, Jac, numIt] = fitModelToData(     ...
    R, model, p0, LB, UB, X, T, data, BCdata, maxIt, plotting, calcJacobian )
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% created 15.02.2017
%
% Notes: -
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for recording the development of chiSqr
P_chiSqr = [];

if size(p0,1) < size(p0, 2)
    p0 = p0';
end

if ~exist('BCdata', 'var')
    BCdata = [];
end

if ~exist('plotting', 'var')
    plotting = 0;
end

%% optimization settings
tolFun = 1e-8;
tolX   = 1e-8;
diffMinChange = 0;
diffMaxChange = Inf;

options = optimoptions(@lsqnonlin, ...
    'Display',            'off', ...  % or 'final'
    'FunctionTolerance',  tolFun,...
    'StepTolerance',      tolX,...
    'DiffMinChange',      diffMinChange,...
    'DiffMaxChange',      diffMaxChange,...
    'MaxIterations',      maxIt, ...
    'Jacobian',           calcJacobian );
%     'FinDiffRelStep',FinDiffRelStep, ...


%% Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pFit, chiSqr, ~, ~, output, ~, Jac] = lsqnonlin(@residual, p0, LB, UB, options);
fprintf('\n')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


numIt = output.iterations;

%% eo fitModelToData



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RESIDUAL function
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [r, J]=residual(p)
        
        amplitude = p(end);
        
        % simulate data
        P = getParameter(model, p, R.lambda);
        
        if nargout == 2
            %% Differential Quotient for Jacobian
            for p_i = 1:np-1 %
                
                p1 = p;
                p2 = p;
                
                dh = max(1e-3*p(p_i), 1e-4);
                p1(p_i) = p(p_i) + dh;
                p2(p_i) = p(p_i) - dh;
                
                P1 = getParameter(model, p1, R.lambda);
                P2 = getParameter(model, p2, R.lambda);
                
                [XC1, TY1, Ysum1] = forwardSimulation(P1, R, model);
                [XC2, TY2, Ysum2] = forwardSimulation(P2, R, model);
                
                [data1] = getFOVdata(XC1, TY1, Ysum1, X, T, p(end));
                [data2] = getFOVdata(XC2, TY2, Ysum2, X, T, p(end));
                
                grad = 1/(2*dh) * (data1 - data2);
                
                J(:, p_i) = grad(:);
            end
            
            %% calc Jacobian entries for amplitude
            [data] = getFOVdata(XC, TY, Ysum, X, T, 1);
            J(:, end) = data(:);
        else
            [XC, TY, Ysum]       = forwardSimulation(P, R, model);
        end
        
        calcData = getFOVdata( XC, TY, Ysum, X, T, amplitude);
        
        % calc residual
        r = calcData(:) - data(:);
        
        if plotting
%             fprintf('.');
            
            subplot(121);
            plot(T, data', 'x')
            completeData = getFOVdata(XC, TY, Ysum, X, TY, amplitude);
            oplot(TY, completeData, '-')
            posMax = gmax(data);
            pName = getParameter(model, 'pInfo');
            pTxt = {numel(p),1};
            for i = 1:length(p)
                pTxt{i} = sprintf('%5s: %2.3f',  pName{i}, p(i));
            end
            text(5, max(posMax*0.8, 0.8), pTxt)
            
            ChiSqr = sum(r(:).^2) / numel(r);
            %             Loss = getLoss_V2(modelNo, p);
            
            theta = p(2) * p(4) / ( p(1)*(p(3) + p(4) ));
            Loss = 1-exp(-theta*10);
            
            pTxt = [];
            pTxt{1} = num2str(ChiSqr, 'ChiSqr: %2.3e');
            pTxt{2} = num2str(Loss, 'Loss: %2.3e');
            text(5, max(posMax*0.3, 0.3), pTxt);
            ylim([-0.1 1.1])
            
            subplot(122)
            P_chiSqr(end+1) = ChiSqr;
            semilogy(P_chiSqr)
            drawnow
        end
        
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EO RESIDUAL function for lsqnonlin
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end