%% Info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains some exemplary code used for studys of experimental
% design for tracer transport in plants. First, some artificial data is
% generated using a forward simulation of the transport model of Buehler et
% al. 2014. This data set is then used for the reconstruction of the 
% original model parameters by inverse modeling. The experimental design
% takes place by defining different windows of the data (simulating 
% interrupted data acquisition). Results are written to a .csv file
%
% Buehler, J., von Lieres, E., and Huber, G. (2014). A class of
% compartmental models for long-distance tracer transport in plants.
% J. Theor. Biol. 341, 131-142. doi: 10.1016/j.jtbi.2013.09.023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
init
h1 = figure(1); clf;

%% Write to file switch
writeToFile = 1;

%% define the model
model = { 'N_3', 'IC_1_gauss', 'u_1', 'a_1_2', 'a_2_1', 'a_2_3', 'x0_1', 'xSigma_1', 'ampl'}; % this corresponds to model M13 from Buehler et al. 2014

%% define dimensions
X = 0:10:100; % mm
dT = 1;
T = 0:dT:180; % min

%% define model parameters
lambda = 0.000567*60; %% decay rate for C-11 in min^-1
v   = 2;       % mm min^-1
e12 = 0.3;     % min^-1
e21 = 0.2;     % min^-1
e23 = 0.05;    % min^-1
x0  = 60;      % mm
sigma = 5;     % mm
amplitude = 1; % a.u.
p = [v e12 e21 e23 x0 sigma amplitude];

%% create prerequisite structure
dx = 1;
R = [];
R.xn  = ceil((X(end)-X(1)) / dx);
R.dx  = 1;
R.solverFun = @ode5;
R.lambda = lambda;
R.Tend = T(end);

%% get parameter structure
P = getParameter(model, p, lambda);

%% create artificial data
[XC, TY, Ysum] = forwardSimulation(P, R, model);
imagesc(TY, XC, Ysum)
xlabel('Time T [min]')
ylabel('Spatial position X [mm]')
title('Tracer distribution')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now use this data set for inverse modeling

% reduce data to spatial positions in X
Ysum = getFOVdata(XC, TY, Ysum, X, TY, amplitude);

% set max time of interest
maxTime = 120;

% set sample handling time
changeTime = 1;

% get boundaries
[~, ~, LB, UB] = getParameter(model, 'pInfo');

% start parameter set
p0 = [1.0 0.1 0.1 0.01 50 10 1];

% fitting settings
plotting = 1;
calcJacobian = 'off';
maxIt = 100;

% counter
counter = 0;
scenarioNumber = 0;

%% create zip file of current code and configuration
if writeToFile
    sessionName = mfilename; %#ok<UNRCH> %%#ok<UNRCH>
    outDir = ['..\..\60 Results\'  sessionName '\'];
    mkdir(outDir)
    zipFilename = [outDir 'code.zip']; % create zip file of code for documentation of results
    createCodeZipFile(zipFilename);
    
    outfilename = [outDir sessionName '.csv' ];
    fid = fopen(outfilename, 'w');
    
    fprintf(fid, 'Scenario number;noise level; windowStart; windowWidth; noOfWindows; distanceBetweenWindows; numel(data); sum(data(:));');
    fprintf(fid, 'sum(relStd 1-4); Mean Number of Samples / h;Cumulated Time of Measurement [min];');
    fprintf(fid, 'numIterations (maxIt=%d);ChiSqr;co-variance matrices;', maxIt);
    fprintf(fid, 'u pFit;mean (n=%d); std; relStd_mean; relStd_std; a12 pFit;mean; std; relStd_mean; relStd_std; a21 pFit;mean; std; relStd_mean; relStd_std; a23 pFit;mean; std; relStd_mean; relStd_std;', NoiseRepNo);
    fprintf(fid, 'x0 pFit;mean; std; relStd_mean; relStd_std; xSigma pFit;mean; std; relStd_mean; relStd_std; amplitude pFit;mean; std; relStd_mean; relStd_std;');
    fprintf(fid, 'Time array T;\n');
end

tic
%% add noise to data
NoiseRepNo = 5;
Data = cell(NoiseRepNo, 1);
rng(123); % set seed for reproducibility
noise_level = 7e-3;
for i = 1:NoiseRepNo
    Data{i} = Ysum + noise_level*randn(size(Ysum));
end

%% loop over discrete widths
for windowWidth =  1:dT:120  % in min
    ww2 = windowWidth + changeTime;

    %% loop over window start position
    for windowStart = 20 : dT : 40
        
        %% loop over number of windows
        NoOfWindows = floor((maxTime-windowStart)/ww2) : -1 : 1;
        for noOfWindows = NoOfWindows
            
            %% loop over varying distances between the windows
            if noOfWindows > 1
                Distances = ww2:ww2:floor((maxTime-windowStart-(noOfWindows-1)*ww2)/(noOfWindows-1));
            else % noOfWindows == 1
                Distances = 5000;
            end
            for distanceBetweenWindows = Distances
                
                % set scenario number
                scenarioNumber = scenarioNumber + 1;

                %% define windows in time
                T = [];
                for i = 1:noOfWindows
                    localWindowStart = windowStart + (i-1) * (ww2 + distanceBetweenWindows);
                    T = [T (localWindowStart) : dT : (localWindowStart + (ww2-changeTime) - dT) ];  % ww2-changeTime = windowWidth
                end
                
                %% create prerequisite structure
                R.Tend = max(gmax(T), 55);
                
                if noise_level > 0
                    PFit = zeros(numel(p0), NoiseRepNo);
                    Loss = zeros(1, NoiseRepNo);
                    RelStdLoss = zeros(1, NoiseRepNo);
                    ChiSqr = zeros(1, NoiseRepNo);
                    RelStd = zeros(numel(p0), NoiseRepNo);
                    MaxEigCov = zeros(1, NoiseRepNo);
                    MinEigCov = zeros(1, NoiseRepNo);
                    NumIt  = zeros(1, NoiseRepNo);
                    Cov_m = cell(NoiseRepNo,1);
                    numelData = zeros(1, NoiseRepNo);
                    sumData = zeros(1, NoiseRepNo);
                else
                    PFit = zeros(numel(p0), 1);
                    ChiSqr = zeros(1, 1);
                    RelStd = zeros(numel(p0), 1);
                    DetCov = zeros(1, 1);
                    MaxEigCov = zeros(1, 1);
                    MinEigCov = zeros(1, 1);
                    NumIt  = zeros(1, 1);
                    Cov_m = cell(1,1);
                end
                
                t = tic;
%                 parfor noiseRepetition = 1:NoiseRepNo % or alternatively non-parallel:
                for noiseRepetition = 1:NoiseRepNo
                    
                    counter = counter + 1;
                    data = getFOVdata(X, TY, Data{noiseRepetition}, X, T, amplitude);
                     
                    plot(T, data, '-x');
                    xlim([0 180]);
                    ylim([-0.1 1.1]);
                    title(['Design ' num2str(scenarioNumber)]);
                    drawnow;

                    %% fit IC data
                    [pFit, chiSqr, Jac, numIt] = fitModelToData(R, model, p0, LB, UB, X, T, data, [], maxIt, plotting, calcJacobian);
                    
                    %% calc statistical measures
                    nmp = numel(data) - numel(pFit);
                    ssq = chiSqr / nmp;
                    cov_m = full(inv(Jac'*Jac)); % co-variance matrix
                    pcov = ssq * cov_m;
                    stdDev_cov = sqrt(diag(pcov));
                    
                    % calculate relative std deviation
                    relStd = stdDev_cov ./ pFit * 100;
                    
                    %% disp basic info
                    fprintf('noise: %f, win width: %d, win start: %d, no win: %d, numIt: %d, fit u: %f\n', noise_level, windowWidth, windowStart, noOfWindows, numIt, pFit(1));
                    
                    numelData(noiseRepetition) = numel(data);
                    sumData(noiseRepetition) = sum(data(:));
                    
                    
                    PFit(:, noiseRepetition)    = pFit;
                    ChiSqr(noiseRepetition)     = chiSqr;
                    RelStd(:, noiseRepetition)  = relStd;
                    NumIt(noiseRepetition)      = numIt;
                    Cov_m{noiseRepetition}      = cov_m;
                end
                toc(t);
                
                pFit_mean = mean(PFit, 2);
                pFit_std  = std(PFit, 0, 2);
                chiSqr_mean = mean(ChiSqr);
                chiSqr_std  = std(ChiSqr);
                relStd_mean = mean(RelStd, 2);
                relStd_std  = std(RelStd, 0, 2);
                
                cumTimeOfMeasurement_min = windowWidth*noOfWindows;
                meanNoSamples_h = 60 / cumTimeOfMeasurement_min;
                                
                %% output pFit, chiSqr, relStd, ...
                if writeToFile

                    fprintf(fid, '%d;%f;%d;%d;%d;%d;%d;%f;%f;%f;', scenarioNumber, noise_level, windowStart, windowWidth, noOfWindows, distanceBetweenWindows, numelData(1), sumData(1));
                    fprintf(fid, '%f;%f;%f;', sum(relStd_mean(1:4)), meanNoSamples_h, cumTimeOfMeasurement_min);
                    
                    fprintf(fid, '%d_', NumIt);
                    fprintf(fid, ';');
                    
                    fprintf(fid, '%f_', ChiSqr);
                    fprintf(fid, ';');
                    
                    for i = 1:numel(Cov_m)
                        fprintf(fid, '%f_', Cov_m{1});
                        fprintf(fid, '#');
                    end
                    fprintf(fid, ';');
                    
                    for i=1:numel(pFit_mean)
                        fprintf(fid, '%f_', PFit(i,:));
                        fprintf(fid, ';');
                        fprintf(fid, '%f;%f;%f;%f;', pFit_mean(i), pFit_std(i), relStd_mean(i), relStd_std(i));
                    end
                    
                    fprintf(fid, '%d;', T);
                    
                    fprintf(fid, '\n');
                    
                end
            end
        end
    end
end

fprintf('Counter: %d\n', counter);
toc

if writeToFile
    fclose(fid);
end

%% Final message
fprintf('---------------------------------------------\n');
fprintf('Succesfully finished\n')
fprintf('---------------------------------------------\n');
%%

