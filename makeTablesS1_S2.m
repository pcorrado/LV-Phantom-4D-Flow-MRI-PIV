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

% This script reproduces supplemental figures S2 and S3 from the
% manuscript, plot bland-altman and correlation plots comparing MRI scans
% (Day 1) vs. the PIV scan (average heartbeat).

clear all; close all; %#ok<CLALL>

% UPDATE THESE TO MATCH YOUR FILE STRUCTUE
BASE_PIV_DIR = "E:\PIV_Data\Registered to Day 1";
BASE_MRI_DIR = "E:\MRI_Data\Day 1";

trajectories = ["PC VIPR", "PC SOS"];
scanTime = ["30","20","10","5","2.5"];
POINTS_TO_PLOT = 20000;

f1 = figure();
f1.Position=[400,-600,550,1300];
f2 = figure();
f2.Position=[400,-600,550,1300];

% Loop through scan times
for row=1:numel(scanTime)
    %Loop through trajectories
    for col=1:numel(trajectories)
        fprintf('Analyzing row %i of %i, column %i of %i.\n',row,numel(scanTime),col,numel(trajectories));
    
        traj = trajectories(col);
        pivDir = sprintf('%s\\%s\\Average Beat',BASE_PIV_DIR,traj);
        maskFile = fullfile(pivDir,'..\MASK.nii');
        mriDir = sprintf('%s\\%s\\%s_Minutes\\',BASE_MRI_DIR,traj,scanTime(row));
        
        % Get average and difference of points in the masks in order to
        % generate the plots
        [r,bias,loa1,loa2,~,~,avg,diff] = computeSpeedBlandAltman(maskFile,pivDir,mriDir,0);
        ind = randperm(numel(avg),POINTS_TO_PLOT);
        avg = avg(ind)/10;
        diff = diff(ind)/10;
        
        % Plot bland-altman graph
        figure(f1);
        subplot('Position',[(col-1)*.495+0.12,(5-row)*0.195+0.05,0.36,0.145]);
        scatter1 = scatter(avg,diff,10,[0.3,0.3,0.3],'o','filled'); 
        scatter1.MarkerFaceAlpha = 0.25;
        scatter1.MarkerEdgeAlpha = 0.25;
        hold on;
        plot([0,15],repmat(bias/10,1,2),'-b');
        plot([0,15],repmat(loa1/10,1,2),'--r');
        plot([0,15],repmat(loa2/10,1,2),'--r');
        hold off;
        ylim([-15,15]);
        xlim([0,15]);
        xticks([0,5,10,15]);
        yticks([-15,-7.5,0,7.5,15]);
        xlabel('(MRI speed + PIV speed)/2 (cm/s)');
        ylabel({'MRI speed - PIV speed', '(cm/s)'});
        
        % Plot correlation graph
        figure(f2);
        subplot('Position',[(col-1)*.5+0.1,(5-row)*0.2+0.04,0.35,0.15]);
        scatter1 = scatter(avg-diff/2,avg+diff/2,10,[0.3,0.3,0.3],'o','filled'); 
        scatter1.MarkerFaceAlpha = 0.25;
        scatter1.MarkerEdgeAlpha = 0.25;
        hold on;
        plot([0,20],[0,20],':k','LineWidth',2);     
        xlabel('MRI speed (cm/s)');
        ylabel('PIV speed (cm/s)');
        p = polyfit(avg-diff/2,avg+diff/2,1);
        rtxt = sprintf('r=%.2f',r);
        eq = sprintf('y = %.2fx + %.2f',p(1),p(2));
        text(0.5,17.5,{rtxt,eq},'FontSize',9);
        plot(0:20,p(1)*(0:20)+p(2),'--b','LineWidth',1);
        hold off;
         ylim([0,20]);
        xlim([0,20]);
        xticks([0,5,10,15,20]);
        yticks([0,5,10,15,20]);
    end
end
