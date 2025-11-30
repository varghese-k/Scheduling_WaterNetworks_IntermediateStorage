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

A=[Q(1:8,1:256)];
C=[Q(11,1:256)];
b=10^6.*[.42; .5133; .667; .867; .4067; 2.053; .6733; .48];

[x,f,itn] = linprog(A,b,C',0,100,10^-10)
X=x./3600;
t=ones(256,1);
T=t'*X;
