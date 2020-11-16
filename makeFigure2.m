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

% This script reproduces figure 2 from the manuscript, which shows
% magnitude and velocity-map images from MR images at various levels of
% retrospective underampling.
clear all; close all; %#ok<CLALL>

% THESE MUST BE UPDATED TO YOUR FILE STRUCTURE
BASE_PIV_DATA = "E:\PIV_Data\Registered to Day 1";
BASE_MRI_DATA = "E:\MRI_Data\Day 1";

trajectories = ["PC VIPR", "PC SOS"];
scanTime = ["30","20","10","5","2.5"];

% Pick a representative slice and time frame
t = 13;
slice = 32;

% Set image size
imgSize = [144,64];
bigSize = imgSize.*[numel(trajectories),numel(scanTime)]+[5,10];
bigImg = ones(bigSize);
yInd = [1,72;1,72];
xInd = [1,64;1,64];

% Loop through trajectories and scan times, load images, and put them in
% the right place
for row = 1:numel(trajectories)
    for col=1:numel(scanTime)
        fprintf('Analyzing column %i of %i.\n',col,numel(scanTime));

        traj = trajectories(row);
        maskFile = string(sprintf('%s\\%s\\MASK.nii',BASE_PIV_DATA,traj));
        mriDir = string(sprintf('%s\\%s\\%s_Minutes\\',BASE_MRI_DATA,traj,scanTime(col)));
        
        mask = niftiread(maskFile);
        mag = niftiread(fullfile(mriDir,'MAG.nii'));
        velZ = niftiread(fullfile(mriDir,'VELZ.nii'));
        mag = fliplr(rot90(squeeze(double(mag(xInd(row,1):xInd(row,2),slice,yInd(row,1):yInd(row,2),t))),1));
        mag = double(mag-0)/30000;
        velZ = fliplr(rot90(squeeze(double(velZ(xInd(row,1):xInd(row,2),slice,yInd(row,1):yInd(row,2),t))),1));
        velZ = double(velZ+100)/200;
        
        bigImg(((row-1)*(imgSize(1)+5)+1):(row*imgSize(1)+5*(row-1)),((col-1)*(imgSize(2)+2)+1):(col*imgSize(2)+2*(col-1))) = [mag;velZ];

    end
end
imshow(bigImg,[0,1]);