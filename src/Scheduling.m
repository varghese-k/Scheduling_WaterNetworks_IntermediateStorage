%--------------------------------------------------------------------------
% Program to create a schedule for the Arcot network
%--------------------------------------------------------------------------
% Clear workspace
home; clear; ...clear import; 
    path(pathdef); rng('default'); format short;
close all;
date_str = string(datetime("now","Format","yyMMdd"));
rightnow = string(datetime("now","Format","HHmm"));


% Set path
[currentFolder,~,~] =... 
    fileparts([mfilename('fullpath') '.m']);
workingDir = fileparts(currentFolder);
srcDir = fullfile(workingDir,'src');
epanetDir = fullfile(workingDir,'epanet');
outDir = fullfile(workingDir,'out');
outputDir = fullfile(outDir,date_str,rightnow);
mkdir(outputDir);
addpath(...
    workingDir,...
    outDir,...
    currentFolder,...
    epanetDir);   %

%--------------------------------------------------------------------------
                % Loading Epanet library and Epanet file %
%--------------------------------------------------------------------------
% load Epanet library 
loadlibrary('epanet2','epanet2.h');
% Load Epanet input file (.inp format)
epanetloadfile('Network_Data.inp');
 
%--------------------------------------------------------------------------
                        % Hydraulic Simulation
%--------------------------------------------------------------------------

calllib('epanet2','ENopenH');

% The following are the index of the pipes leading to the inlet of the
% eight tanks in the Arcot network. Please note, these are not pipe numbers
pipeIndices = [7;8;10;11;18;17;15;13];

Q = zeros(11,256);

% Binary variables to set the status of pipes (ON/OFF)
j1 = 0;
j2 = 0;
j3 = 0;
j4 = 0;
j5 = 0;
j6 = 0;
j7 = 0;
j8 = 0;
j9 = 0;
j = [j1; j2;j3;j4;j5;j6;j7;j8];
m = [128;64;32;16;8;4;2;1];
nTanks = length(m);
for j1=0:1
    for j2=0:1
        for j3=0:1
            for j4=0:1
                for j5=0:1
                    for j6=0:1
                        for j7=0:1
                            for j8=0:1
                                err = calllib('epanet2','ENinitH',0);
                                j = [j1; j2;j3;j4;j5;j6;j7;j8];
                                for i = 1:8
                                    err = calllib('epanet2','ENsetlinkvalue',pipeIndices(i),11,j(i));
                                end
                                [err,time] = calllib('epanet2','ENrunH',0);
                                for i = 1:8
                                    [err,Q(i,m'*j+1)] = calllib('epanet2','ENgetlinkvalue',pipeIndices(i),8,Q(i,m'*j+1));
                                end
                                [err,Q(9,m'*j+1)] = calllib('epanet2','ENgetlinkvalue',20,13,Q(9,m'*j+1));
                                [err,Q(10,m'*j+1)] = calllib('epanet2','ENgetlinkvalue',20,8,Q(10,m'*j+1));
                                Q(11,m'*j+1) = Q(9,m'*j+1)*100/(-0.0069*Q(10,m'*j+1)^2+1.5972*Q(10,m'*j+1)-11.667);
                            end
                        end
                    end
                end
            end
        end
    end
end
calllib('epanet2','ENclose');
unloadlibrary('epanet2');

%--------------------------------------------------------------------------
                        % Linear Program to find optimal states
%--------------------------------------------------------------------------

Aeq = [Q(1:nTanks,1:256)]*3600; % Convert flow rates from L/s to L/h
C = [Q(nTanks+3,1:256)];
C(1) = 0.1; % Small positive cost for when all valves are off 
beq = 10^6.*[.63; .77; 1; 1.3; .61; 3.08; 1.01; .72]; % Cumulative demand in L

[x,f,exitFlag] = linprog(C',[],[],[Aeq;ones(1,256)],[beq;24],zeros(size(Aeq,2),1),24*ones(size(Aeq,2),1));

t = ones(256,1);
T = t'*x; % Total pumping time


%--------------------------------------------------------------------------
                        % Integer program to find schedule
%--------------------------------------------------------------------------


prob = optimproblem; 

% Withdrawal pattern for all tanks (fraction of average withdrawal rate in each hour
% starting at midnight)
hours_demand = 1:1:24;
W = [0.2, 0.2, 0.2, 0.2, 0.4, 0.8, 2.25, 2.25, 2.25, 2.25, 1, 1, 0.6, 2.25, 0.7, 0.7, 0.7, 2.25, 0.9, 0.9, 0.7, 0.7, 0.4, 0.2];
W_cumulative = cumsum(W)/length(W); % Fraction of the daily demand completed at the end of each hour, assumed to be same across the thanks

Capacity = beq/3;  % Capacity reduced to one-third of the daily demand  based on design requirements. Actual capacities were larger.

% Active states and times
A = Aeq(:,x>0);
t = x(x>0);

numStates = numel(t);
minSliceDuration = 0.5; % min duration of each slice in hours

% Slice the total duration of each state into smaller segments

nSlices = zeros(numStates,1);
tau     = zeros(numStates,1);  % per-slice duration (hours) for each state

for m = 1:numStates
    % maximum integer number of equal slices with each >= 0.5 h
    nSlices(m) = max(1, floor(t(m) / minSliceDuration));

    tau(m) = t(m) / nSlices(m);
end

% total number of slices across all states
totalSlices = sum(nSlices); % This + 1 give  the number of event points in the day. 

% Define integer program 

% Volume of water in each tank at each time point (upper and lower bounds
% are 10% and 90% of tank capacity)
V = optimvar("V",nTanks, totalSlices,Type="continuous",LowerBound=repmat(Capacity,1,totalSlices)*0.05,UpperBound=repmat(Capacity,1,totalSlices)*0.95);

% Binary variable indicating whether state i is active in time interval j
y = optimvar("y",numStates, totalSlices,Type="integer",LowerBound=0,UpperBound=1);

% Continuous variable describing the end points of each time slot in the
% schedule
Tvals = optimvar("Tvals",1, totalSlices, Type="continuous",LowerBound=0,UpperBound=sum(t)+0.1);

% Cumulative demand at each event point (given as a fraction of the total
% daily demand)
D = optimvar('D', 1,  totalSlices, Type="continuous", LowerBound=0);  % interpolation weights

% variable for interpolating the hourly demand vector to find teh
% cumulative demand at each event point
lambda = optimvar('lambda', 24,  totalSlices, 'LowerBound', 0);  % interpolation weights

%%%%% The follwoing lines describes the constraints to define a
%%%%% piecewise-linear interpolation. Required because the event points in
%%%%% the schedule do not necessarily coincide with the hours of the day,
%%%%% for which the demands are availble. 

% Binary selector z(i,k): interval m active for time slot j 
z = optimvar('z', length(hours_demand), totalSlices, 'Type','integer','LowerBound',0,'UpperBound',1);

% Each slice must choose exactly one interval
prob.Constraints.oneInterval = sum(z,1) == 1;

% Big-M constants
M_t  = max(hours_demand) - min(hours_demand) + 10;  % safe big-M for time
M_W  = max(W_cumulative) - min(W_cumulative) + 10;  % safe big-M for W


intervalConstr = optimconstr(4*(length(hours_demand))*totalSlices);
row = 1;
for i = 0:length(hours_demand)-1
    if i == 0
        x0 = 0;
        y0 = 0;
    else
        x0 = hours_demand(i);
        y0 = W_cumulative(i);
    end
    x1 = hours_demand(i+1); 
    y1 = W_cumulative(i+1);

    slope = (y1 - y0) / (x1 - x0);

    for j = 1:totalSlices

        % Bound t_val inside interval
        intervalConstr(row)   = Tvals(j) >= x0 - M_t*(1 - z(i+1,j)); row = row+1;
        intervalConstr(row)   = x1 + M_t*(1 - z(i+1,j)) >= Tvals(j); row = row+1;

        % Interpolated value (lower and upper)
        intervalConstr(row)   = D(j) >= y0 + slope*(Tvals(j) - x0) - M_W*(1 - z(i+1,j)); row = row+1;
        intervalConstr(row)   = y0 + slope*(Tvals(j) - x0) + M_W*(1 - z(i+1,j)) >= D(j); row = row+1;

    end
end

prob.Constraints.interval = intervalConstr;


for j = 1:totalSlices
    prob.Constraints.("single_state" + j) = sum(y(:,j)) == 1;
    
    if j==1
        prob.Constraints.("Tval_definition" + j) = tau'*y(:,j) == Tvals(j);
    else
        prob.Constraints.("Tval_definition" + j) = tau'*y(:,j) == Tvals(j) - Tvals(j-1);
    end
    
    % cumulative inflow - cumulative ouflow + initial volume == current
    % volume  (as this is a cyclic schedule, initial volume is same as final)
    for k = 1:nTanks
        prob.Constraints.("Volume" + (j-1)*nTanks + k) = (A(k,:).*tau')*y(:,1:j) - D(j)*beq(k) + V(k,totalSlices) == V(k,j) ;
    end
end

for m = 1:numStates
    prob.Constraints.("States_active_time" + m) = sum(y(m,:)) == nSlices(m);
end

sol = prob.solve();