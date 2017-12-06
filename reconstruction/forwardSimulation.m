function [XC, TY, Ysum, Y] = forwardSimulation(P, R, model)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% created 14.02.2017 
% example for R:
% R.xn = ceil(L / dx);
% R.dx = dx;
% R.Tend = Tend;
% R.dm = 'lin3';
% R.rt = 1e-6;
% R.at = 1e-8;
% R.solverFun = @ode15s; % or
% R.solverFun = @ode5;
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get parameter structure
N  = P.N;
%%

dx   = R.dx;
Tend = R.Tend;

xa =  - ceil(gmax(P.x0)  +  10 * gmax(P.xSigma));
xe = R.xn * dx;

XC = (xa : dx : ( xe+dx))';
xn = length(XC);


%% calculate initial function y0
y_0 = initial();

%% exchange matrix
A = P.A';

% matrices of u and d
U = repmat(P.u,xn,1);
D = repmat(P.d,xn,1);

if sum(P.d(:)) > 0
    flag_D = true;
else
    flag_D = false;
end

%% calculate time vector
CFL = .9;
nt = (Tend - 0) * ceil(gmax(P.u)/CFL/dx) + 1;
TY = linspace(0, Tend, nt);


%% CALCULATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y_ = ode5(@right, TY, y_0);

Y       = permute(reshape(Y_(1:nt*xn*N), nt, xn, N), [2 1 3]);
Ysum    = sum(Y,3); % global tracer distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% error checking
if find(isnan(Y))
    error('Solutions contains NaNs')
end
if gmin(Y) < -0.01
    warning(  'forwardSimulation:negativeSolution', 'Solution is negativ!');
end

%% eo forwardSimulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initial functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function y0 = initial()
        y0 = zeros(xn, N);
        for k = 1:N
            switch P.initialFun{k}
                case 'gauss'
                    if (P.IC0(k) > 0)
                        y0(:,k) = 1/(sqrt(2*pi)*P.xSigma(k)) * exp(-((XC+P.x0(k)).^2/(2*P.xSigma(k)^2))) * P.IC0(k);
                    end
                case 'smoothEdge'
                    y0(:,k) = smoothEdge(XC, P.x1(k), P.x2(k), P.edgeSteepNess(k));
                case 'zero'
                    % do nothing
                otherwise
                    error('undefined initial function');
            end
        end
        y0 = y0(:);
    end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% right hand side
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function y_new = right(~,y)
        
        %% forward simulation
        rho  = reshape(y(1:xn*N),xn,N);
        
        % lin3, linear upwind methods of 5th order, periodic
        yd = [rho(end-2,:); rho(end-1,:); rho(end,:); rho; rho(1,:); rho(2,:)];
        
        drhodx =  ( - 1/30*yd(1:xn,:)  + 1/4*yd(2:xn+1,:) - yd(3:xn+2,:) ...
            +  1/3*yd(4:xn+3,:) + 1/2*yd(5:xn+4,:) - 1/20*yd(6:xn+5,:) ) / dx;        
        
        dyc = - U   .* drhodx;
        
        % Diffusion
        if flag_D
            % lin3 periodic
            drhodx2 = (yd(3:xn+2,:)-2*yd(4:xn+3,:)+yd(5:xn+4))/dx^2;
            dyc = dyc + D .* drhodx2;
        end
        
        % Kopplungsterme
        dyc = dyc + rho * A;
        
        y_new = dyc(:);
    end
%% eo right hand side
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end


