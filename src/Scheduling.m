%--------------------------------------------------------------------------
% Program to create a schedule for the water networks with intermediate
% storage
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

%% ------------------------------------------------------------------------
                        % Linear Program to find optimal states
%  ------------------------------------------------------------------------

Aeq = [Q(1:nTanks,1:256)]*3600; % Convert flow rates from L/s to L/h
C = [Q(nTanks+3,1:256)];
C(1) = 0.1; % Small positive cost for when all valves are off 
beq = 10^6.*[.63; .77; 1; 1.3; .61; 3.08; 1.01; .72]; % Cumulative demand in L

[x,f,exitFlag] = linprog(C',[],[],[Aeq;ones(1,256)],[beq;24],zeros(size(Aeq,2),1),24*ones(size(Aeq,2),1));

% t = ones(256,1);
% T = t'*x; % Total pumping time


%% ------------------------------------------------------------------------
                        % Integer program to find schedule
%  ------------------------------------------------------------------------


prob = optimproblem; 

% Withdrawal pattern for all tanks (fraction of average withdrawal rate in each hour
% starting at midnight)
hours_demand = 1:1:24;
W = [0.2, 0.2, 0.2, 0.2, 0.4, 0.8, 2.25, 2.25, 2.25, 2.25, 1, 1, 0.6, 2.25, 0.7, 0.7, 0.7, 2.25, 0.9, 0.9, 0.7, 0.7, 0.4, 0.2];
W_cumulative = cumsum(W)/length(W); % Fraction of the daily demand completed at the end of each hour, assumed to be same across the thanks

Capacity = beq/3;  % Capacity reduced to one-third of the daily demand  based on design requirements. Actual capacities were larger.

% Active states and times
A = Aeq(:,x>0);
t_array = x(x>0);

numStates = numel(t_array);
minSliceDuration = 0.5; % min duration of each slice in hours

% Slice the total duration of each state into smaller segments

nSlices = zeros(numStates,1);
tau     = zeros(numStates,1);  % per-slice duration (hours) for each state

for m = 1:numStates
    % maximum integer number of equal slices with each >= 0.5 h
    nSlices(m) = max(1, floor(t_array(m) / minSliceDuration));

    tau(m) = t_array(m) / nSlices(m);
end

% total number of slices across all states
totalSlices = sum(nSlices); % This + 1 give  the number of event points in the day. 

% Define integer program 

% Volume of water in each tank at each time point (upper and lower bounds
% are 10% and 90% of tank capacity)
V = optimvar("V",nTanks, totalSlices,Type="continuous",LowerBound=repmat(Capacity,1,totalSlices)*0.1,UpperBound=repmat(Capacity,1,totalSlices)*0.9);

% Binary variable indicating whether state i is active in time interval j
y = optimvar("y",numStates, totalSlices,Type="integer",LowerBound=0,UpperBound=1);

% Continuous variable describing the end points of each time slot in the
% schedule
Tvals = optimvar("Tvals",1, totalSlices, Type="continuous",LowerBound=0,UpperBound=sum(t_array));

% Cumulative demand at each event point (given as a fraction of the total
% daily demand)
D = optimvar('D', 1,  totalSlices, Type="continuous", LowerBound=0);  % interpolation weights

%%%%% The follwoing lines describes the constraints to define a
%%%%% piecewise-linear interpolation--required because the event points in
%%%%% the schedule do not necessarily coincide with the hours of the day,
%%%%% for which the demands are availble. 

% Binary selector z(i,k): interval m active for time slot j 
z = optimvar('z', length(hours_demand), totalSlices, 'Type','integer','LowerBound',0,'UpperBound',1);

% Each slice must choose exactly one interval
prob.Constraints.oneInterval = sum(z,1) == 1;

% Big-M constants
M_t  = max(hours_demand) - min(hours_demand) + 100;  % safe big-M for time
M_W  = max(W_cumulative) - min(W_cumulative) + 100;  % safe big-M for W


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

row = 1;
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
        prob.Constraints.("Volume" + row) = sum((A(k,:).*tau')*y(:,1:j)) - D(j)*beq(k) + V(k,totalSlices) == V(k,j) ; row = row+1;
    end
end

for m = 1:numStates
    prob.Constraints.("States_active_time" + m) = sum(y(m,:)) == nSlices(m);
end

[sol,fval,exitflag] = prob.solve();

if exitflag<=0
    return;
end

%% PLot tank volumes

% Extract and plot the solution results
V_solution = sol.V;
D_solution = sol.D;
Tvals_solution = sol.Tvals;
yvals_solution = sol.y;

% Plot of tank volumes at different time intervals
% Number of rows = nTanks/2 (if nTanks is even)
nRows = ceil(nTanks/2);   % safer in case nTanks is odd
nCols = 2;

f = figure(1);

for i = 1:nTanks
    subplot(nRows,nCols,i);
    
    % Plot volume vs. time for tank i [the volume at the end is made the volume at the start as well]
    plot([0,Tvals_solution], [V_solution(i,end), V_solution(i,:)], 'LineWidth', 1.5);
    grid on;
    xticks(0:4:24);
    xlim([0,24]);
    % Axis labeling
    if i > nTanks-nCols
        xlabel('Time (h)');
    end
    ylabel('Volume (L)');
    title(sprintf('Tank %d', i));
    
    % Apply y-axis limits
    ylim([0, Capacity(i)]);
end

sgtitle('Tank volumes over the scheduling horizon (y-axis limits are the capacities)');

saveas(f,fullfile(outputDir,'TankVolumes.fig'));

%% Print schedule as a chart

% Active state per interval
[~, activeState] = max(yvals_solution, [], 1);

% Column widths
dt = zeros(1, totalSlices);
dt(1) = Tvals_solution(1);
for j = 2:totalSlices
    dt(j) = Tvals_solution(j) - Tvals_solution(j-1);
end

% X positions
x_left  = [0, cumsum(dt(1:end-1))];
x_right = cumsum(dt);
totalWidth = sum(dt);

% Reverse tank order (Tank 1 at TOP)
tankYpos = nTanks:-1:1;   % row 8 bottom → row 1 top

f = figure(2); hold on;

% Colors for states
cmap = lines(numStates);
cmap = cmap+(1-cmap)*0.5;
% -----------------------------------------------------------------------
% 1) DRAW COLORED CELLS INSIDE GRID FOR EACH TANK
% -----------------------------------------------------------------------
for j = 1:totalSlices
    s = activeState(j);
    col = cmap(s,:);

    for t = 1:nTanks
        if A(t, s) ~= 0
            % y-position for Tank t
            row = tankYpos(t) - 1;
            rectangle('Position', [x_left(j), row, dt(j), 1], ...
                'FaceColor', col,'EdgeColor','none');
        end
    end
   
end

    

% -----------------------------------------------------------------------
% 2) ADD STATE LABELS INSIDE THEIR OWN TOP ROW CELL
% -----------------------------------------------------------------------
stateRowY = nTanks;   % state labels are ABOVE tank 1 row

for j = 1:totalSlices
    s = activeState(j);
    % Draw the cell for the state label
    rectangle('Position', [x_left(j), stateRowY, dt(j), 1], ...
        'FaceColor',[0.95 0.95 0.95],  'EdgeColor','none');
    
    % Text inside the cell
    text(x_left(j) + dt(j)/2, stateRowY + 0.5, sprintf('S%d', s), ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'FontWeight','bold');
end



% -----------------------------------------------------------------------
% 3) ADD TANK LABELS INSIDE FIRST COLUMN
% -----------------------------------------------------------------------
firstColWidth = dt(1);

% Add color for the top left cell
rectangle('Position', [-2.5*firstColWidth, stateRowY, 2.5*firstColWidth, 1], ...
    'FaceColor',[0.95 0.95 0.95],  'EdgeColor','none');

for t = 1:nTanks
    row = tankYpos(t) - 1;

    rectangle('Position', [-2.5*firstColWidth, row, 2.5*firstColWidth, 1], ...
        'FaceColor', [0.95, 0.95, 0.95], 'EdgeColor','none');
    text(-2*firstColWidth, row+0.5, sprintf('Tank%d', t), ...
        'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
        'FontWeight','bold');
end
     
% -----------------------------------------------------------------------
% 4) ADD GRID LINES AROUND *ALL* CELLS (INCLUDING STATE ROW)
% -----------------------------------------------------------------------
% Horizontal grid lines
for r = 0:(nTanks+1)   % +1 for state row
    line([-2.5*firstColWidth totalWidth], [r r], 'Color','k', 'LineWidth',1.2);
end


% Vertical grid lines
for j = 1:totalSlices
    line([x_left(j) x_left(j)], [0 nTanks+1], 'Color','k', 'LineWidth',1.2);
end
line([totalWidth totalWidth], [0 nTanks+1], 'Color','k', 'LineWidth',1.2);

line([-2.5*firstColWidth -2.5*firstColWidth], [0 nTanks+1], 'Color','k', 'LineWidth',1.2);

% -----------------------------------------------------------------------
% Axis adjustments
% -----------------------------------------------------------------------
xlim([-2.5*firstColWidth totalWidth]);
ylim([0 nTanks+1]);

% X-axis ticks every 4 hours from 0 to 24
xticks(0:4:24);
xticklabels(string(0:4:24));

yticks([]);   % labels drawn manually

xlabel('Time (hours)');
title('Gantt chart showing the schedule');

hold off

saveas(f,fullfile(outputDir,'GanttChart.fig'));


%% Write results to an xlsx file

% Statedata
tankLabels = strcat("Tank ", string(1:nTanks), " (L/h)" )';
tankLabels = [tankLabels; "Total time (h)"];

stateNames = strcat("S", string(1:numStates));

SelectedStates = array2table([A; t_array'], 'VariableNames', stateNames, 'RowNames', tankLabels);

writetable(SelectedStates, fullfile(outputDir, 'Solution.xlsx'), 'Sheet', 'StateData', 'WriteRowNames', true);

% Schedule data (volume)
rowLabels = strcat("Tank ", string(1:nTanks), " (kL)" )';

% Active states
stateNames = strcat("S", string(activeState));

endTimes = strings(size(Tvals_solution));
for k = 1:numel(endTimes)
    endTimes(k) = sprintf('%.2f', Tvals_solution(k));   
end

tableTitle = [strcat("State: ", stateNames) + newline + strcat("End time (h): ", endTimes)];

% Create table for schedule data

scheduleData = array2table(V_solution/1000, 'VariableNames', tableTitle, 'RowNames', rowLabels);


writetable(scheduleData, fullfile(outputDir, 'Solution.xlsx'), 'Sheet', 'ScheduleData', 'WriteRowNames', true);