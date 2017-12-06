function init()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% created 14.02.2017 by Jonas Buehler, Gregor Huber
%
% call this function once to initialize your working environment
%  - add paths of project structure to matlabs path variable
%  - initialize parallel computing environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% add paths from local code structure
fid = fopen('init_paths.txt', 'r');
tline = fgetl(fid);
while ischar(tline)
    addpath(tline);
    tline = fgetl(fid);
end
fclose(fid);

return

%% set some PDE solver parameter
dx = 1;
dt = 1;
dm = 'lin3'; % 

solverFun = @ode15s;



