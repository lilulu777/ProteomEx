% Measurement Error Compiling

%% ===== Data Compiling =====

% Use the first 3 sections to compile RMS error data from separate samples
% Then, run the following code to get RMS curve with SD shading.

%% Initialization: Run this section once
clear;
measData = struct([]);

%% Compiling: Run this section N times, for each of the N samples

% Step 1: Load saved workspace from MasterScript
% Step 2: Edit the following parameters and run this section
n = 1
data_name = 'Workspace_After_ErrorAnalysis_MT1.mat'
load(data_name)

measData(n).Name = data_name;
measData(n).measResults = imgPair.measResults;
measData(n).pixelWidth = imgPair.info.pixelWidth;

n = 2
data_name = 'Workspace_After_ErrorAnalysis_MT2.mat'
load(data_name)

measData(n).Name = data_name;
measData(n).measResults = imgPair.measResults;
measData(n).pixelWidth = imgPair.info.pixelWidth;

n = 3
data_name = 'Workspace_After_ErrorAnalysis_MT3.mat'
load(data_name)

measData(n).Name = data_name;
measData(n).measResults = imgPair.measResults;
measData(n).pixelWidth = imgPair.info.pixelWidth;



%% Saving: Run this section once
clearvars -except measData
save measData_compiled.mat;

%% ===== Data Visualization =====
clear; clc; close all;
load measData_compiled.mat;

%% individual plots
figure(1); clf;

panel_config = [2 3];

movavgWindow = 10;
do_movavg = 0;

axis_limits = [0 1500 0 200];

for i = 1:length(measData)
    
    % Basic calc from original script
    measResults = measData(i).measResults;
    xmin = 1;%in pixels
    xmax = size(measResults,1);
    xvals = [xmin:xmax];
    notNaNindices = xvals(~isnan(measResults([xmin:xmax],3)));
    pixelWidth_pre = measData(i).pixelWidth;
    xvals = notNaNindices * pixelWidth_pre;
    
    measLength = xvals;
    RMSerror = measResults(notNaNindices,3)*pixelWidth_pre;
    
    if do_movavg
        RMSerror = movmean(RMSerror,movavgWindow);
    end
    
    subplot(panel_config(1),panel_config(2),i);
    plot( measLength, RMSerror, 'LineWidth', 2)
    
    axis(axis_limits);
end

%% AVG and STD plot

num_xval = 200;
measResults_all = zeros(length(measData),num_xval);

for i = 1:length(measData)
    
    measResults = measData(i).measResults;
    xmin = 1;%in pixels
    xmax = size(measResults,1);
    xvals = [xmin:xmax];
    notNaNindices = xvals(~isnan(measResults([xmin:xmax],3)));
    pixelWidth_pre = measData(i).pixelWidth;
    xvals = notNaNindices * pixelWidth_pre;
    
    measLength = xvals(1:num_xval);
    RMSerror = measResults(notNaNindices,3)*pixelWidth_pre;
    
    if do_movavg
        RMSerror = movmean(RMSerror,movavgWindow);
    end
    
    measResults_all(i,:) = RMSerror(1:num_xval);
end

measResults_AVG = mean(measResults_all,1);
measResults_STD = std(measResults_all,0,1);

figure(2); clf;
errorbar(measLength, measResults_AVG, measResults_STD);
axis(axis_limits);

figure(3);clf;
h = shadedErrorBar(measLength, measResults_AVG, measResults_STD,'lineprops','b');
set(h.mainLine,'LineWidth',2)
axis(axis_limits);
xlabel('Measurement Length (\mum)','fontsize',12)
ylabel('RMS Error (\mum)','fontsize',12)
    