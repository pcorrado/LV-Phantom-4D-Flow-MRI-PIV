% Copyright (c) 2020 Philip A Corrado, University of Wisconsin-Madison
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% This script reproduces figure 4 from the manuscript, which plots velocity
% RMSE for MRI images relative to PIV images for all scan times and
% trajectories.

clear all; close all; %#ok<CLALL>

% UPDATE THESE TO MATCH YOUR FILE STRUCTURE
BASE_PIV_DATA = "E:\PIV_Data";
BASE_MRI_DATA = "E:\MRI_Data";

piv_days = ["Registered to Day 1"; "Registered to Day 2"];
trajectories = ["PC VIPR", "PC SOS"];

scanTime = ["30","20","10","5","2.5"];
mri_days = ["Day 1", "Day 2"];

% Loop through trajectories and scan times and compute RMSE
for row=1:(numel(trajectories)*numel(scanTime))
    traj = trajectories(ceil(row/numel(scanTime)));
    for day = 1:numel(piv_days)
        fprintf('Computing traj=%s, scanTime=%s, %s.\n', traj,scanTime(mod(row-1,numel(scanTime))+1),mri_days(day));
        PIV_Dirs(row,day) = string(sprintf('%s\\%s\\%s\\%s\\',BASE_PIV_DATA,piv_days(day),traj,"Average Beat")); %#ok<SAGROW>
        MRI_Dirs(row,day) = string(sprintf('%s\\%s\\%s\\%s_Minutes\\',BASE_MRI_DATA,mri_days(day),traj,scanTime(mod(row-1,numel(scanTime))+1))); %#ok<SAGROW>
        vRMSE(row,day) = compareFlowDataSets(fullfile(PIV_Dirs(row,day),'..\MASK.nii'),PIV_Dirs(row,day),MRI_Dirs(row,day),0); %#ok<SAGROW> 
    end
end
vRMSE = vRMSE/10; % Convert from mm/s to cm/s

% Plot result
time = str2double(scanTime);
f = figure(1);
hold on;
plot(time,vRMSE(1:5,1),'-ro','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1,'DisplayName', 'PC VIPR Day 1');
plot(time,vRMSE(1:5,2),'--ro','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1,'DisplayName', 'PC VIPR Day 2');
plot(time,vRMSE(6:10,1),'-ko','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1,'DisplayName', 'PC SOS Day 1');
plot(time,vRMSE(6:10,2),'--ko','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1,'DisplayName', 'PC SOS Day 2');
hold off;
set(gca,'fontsize',13);
t=title('RMS Velocity Error Relative to PIV');
t.FontSize = 17;
xl=xlabel('Scan time (min.)');
yl=ylabel('RMS velocity Error (cm/s)');
xl.FontSize=15;
yl.FontSize=15;
l = legend();
l.String = l.String(1:4);
l.FontSize=14;
f.Position = [440   106   643   692];

