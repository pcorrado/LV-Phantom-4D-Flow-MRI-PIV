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

% This script reproduces figure 3 from the manuscript, which shows
% velocity-map images from PIV and selected MR images for a characteristic
% slice.

clear all; close all; %#ok<CLALL>

% UPDATE THESE TO MATCH YOUR FILE STRUCTURE
imgDirs = ["E:\PIV_Data\Registered to Day 1\PC VIPR\Average Beat";...
           "E:\MRI_Data\Day 1\PC VIPR\20_Minutes";...
           "E:\MRI_Data\Day 1\PC SOS\20_Minutes";...
           "E:\MRI_Data\Day 1\PC VIPR\5_Minutes";...
           "E:\MRI_Data\Day 1\PC SOS\5_Minutes"];
       
maskFiles = ["E:\PIV_Data\Registered to Day 1\PC VIPR\MASK.nii";...
             "E:\PIV_Data\Registered to Day 1\PC VIPR\MASK.nii";...
             "E:\PIV_Data\Registered to Day 1\PC SOS\MASK.nii";...
             "E:\PIV_Data\Registered to Day 1\PC VIPR\MASK.nii";...
             "E:\PIV_Data\Registered to Day 1\PC SOS\MASK.nii"];

t = 13;
slice = 34;

imgSize = [96,27];
bigSize = imgSize.*[1,5];
bigImg = ones(bigSize);
 
% Loop through images, crop image, trace boundary, and add to big figure.
for col=1:numel(imgDirs)
    fprintf('Analyzing column %i of %i.\n',col,numel(imgDirs));

    mask = niftiread(maskFiles(col));

    velX = niftiread(fullfile(imgDirs(col),'VELX.nii')).*mask;
    velY = niftiread(fullfile(imgDirs(col),'VELY.nii')).*mask;
    velZ = niftiread(fullfile(imgDirs(col),'VELZ.nii')).*mask;

    velX = fliplr(rot90(squeeze(double(velX(12:38,slice,14:45,t))),1));
    velX = double(velX+150)/300;
    velY = fliplr(rot90(squeeze(double(velY(12:38,slice,14:45,t))),1));
    velY = double(velY+150)/300;
    velZ = fliplr(rot90(squeeze(double(velZ(12:38,slice,14:45,t))),1));
    velZ = double(velZ+150)/300;
    
    mask = fliplr(rot90(squeeze(double(mask(12:38,slice,14:45,t))),1));
    
    BD = bwboundaries(mask);
    contours{col} = BD{1}; %#ok<SAGROW>
    bigImg(:,((col-1)*(imgSize(2))+1):(col*imgSize(2))) = [velX;velZ;velY];
end
imshow(bigImg,[0,1]);
hold on;
for col = 1:numel(contours)
    plot(contours{col}(:,2)+(col-1)*(imgSize(2)), contours{col}(:,1), '--r', 'LineWidth', 1.5);
    plot(contours{col}(:,2)+(col-1)*(imgSize(2)), contours{col}(:,1)+(imgSize(1))/3, '--r', 'LineWidth', 1.5);
    plot(contours{col}(:,2)+(col-1)*(imgSize(2)), contours{col}(:,1)+((imgSize(1))*2)/3, '--r', 'LineWidth', 1.5);
end
hold off;
