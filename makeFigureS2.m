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

% This script reproduces supplemental figure S4 from the manuscript, which
% tabulates RMSE for MR images for each PIV data set available: each of the
% 4 hearbeats recorded and the average heartbeat.

BASE_PIV_DIR = "E:\PIV_Data";
BASE_MRI_DIR = "E:\MRI_Data";

piv_days = ["Registered to Day 1"; "Registered to Day 2"];
beats = ["Heartbeat 1", "Heartbeat 2", "Heartbeat 3", "Heartbeat 4", "Average Beat"];

trajectories = ["PC VIPR", "PC SOS"];
scanTime = ["30","20","10","5","2.5"];
mri_days = ["Day 1", "Day 2"];

T = table('Size',[10,7],'VariableTypes',{'string','string','string','string','string','string','string'},...
'VariableNames',{'Scan','Scan_Time','Beat_1','Beat_2','Beat_3','Beat_4','Average_Beat'});

% Loop through PIV and MRI datasets and fil in RMSE table
for col=1:numel(beats)
    for row=1:(numel(trajectories)*numel(scanTime))
        traj = trajectories(ceil(row/numel(scanTime)));
        T(row,1) = {traj};
        T(row,2) = {scanTime(mod(row-1,numel(scanTime))+1)};
        for day = 1:numel(piv_days)
            PIV_Dirs(row,col,day) = string(sprintf('%s\\%s\\%s\\%s\\',BASE_PIV_DIR,piv_days(day),traj,beats(col))); %#ok<SAGROW>
            MRI_Dirs(row,col,day) = string(sprintf('%s\\%s\\%s\\%s_Minutes\\',BASE_MRI_DIR,mri_days(day),traj,scanTime(mod(row-1,numel(scanTime))+1))); %#ok<SAGROW>
            vRMSE(row,col,day) = compareFlowDataSets(fullfile(PIV_Dirs(row,col,day),'..\MASK.nii'),PIV_Dirs(row,col,day),MRI_Dirs(row,col,day),0); %#ok<SAGROW> 
        end
        T(row,col+2) = {sprintf('%2.1f (%2.1f)',mean(vRMSE(row,col,:),3)/10,std(vRMSE(row,col,:),[],3)/10)}; % Divide by 10 to convert from mm/s to cm/s
    end
end
