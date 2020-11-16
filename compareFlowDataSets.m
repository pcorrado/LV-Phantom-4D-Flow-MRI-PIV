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

function vRMSE = compareFlowDataSets(maskNiiFile,referenceImgNiiDir,testImgNiiDir,tShift)
% This function compares two time-resolved, volumetric velocity datasets
% (i.e. PIV vs. MRI or MRI vs. MRI) voxel-by-voxel, and returns the
% root-mean-square error in the velocity vector of the test image relative
% to the reference image. The output vRMSE is in mm/s.

    % Load images
    maskImg = niftiread(maskNiiFile);
    maskImg = logical(maskImg);
    refVelX = niftiread(fullfile(referenceImgNiiDir,'VELX.nii'));
    refVelY = niftiread(fullfile(referenceImgNiiDir,'VELY.nii'));
    refVelZ = niftiread(fullfile(referenceImgNiiDir,'VELZ.nii'));
    testVelX = circshift(niftiread(fullfile(testImgNiiDir,'VELX.nii')),[0,0,0,tShift]);
    testVelY = circshift(niftiread(fullfile(testImgNiiDir,'VELY.nii')),[0,0,0,tShift]);
    testVelZ = circshift(niftiread(fullfile(testImgNiiDir,'VELZ.nii')),[0,0,0,tShift]);
    
    if (size(refVelX,4)~=size(testVelX,4))
        error('Reference and test images have different number of time frames.');
    end
    
    dVelSquaredTotal = []; % container to hold the square velocity error for all voxels and time frames
    % Loop through time frames
    for tt = 1:size(refVelX,4)
       diffX = testVelX(:,:,:,tt)-refVelX(:,:,:,tt); 
       diffY = testVelY(:,:,:,tt)-refVelY(:,:,:,tt); 
       diffZ = testVelZ(:,:,:,tt)-refVelZ(:,:,:,tt); 
       dVelSquared = double(diffX).^2 + double(diffY).^2 + double(diffZ).^2;
       dVelSquaredTotal = [dVelSquaredTotal;dVelSquared(maskImg(:,:,:,tt)==1)]; %#ok<AGROW>
    end
    vRMSE = sqrt(mean(dVelSquaredTotal)); % Take the root and the mean of the square
end

