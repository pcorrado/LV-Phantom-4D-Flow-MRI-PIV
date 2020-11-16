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
function [r,bias,loa1,loa2,CoV,p,avg,diffRT] = computeKEBlandAltman(maskNiiFile,referenceImgNiiDir,testImgNiiDir,tShift)
    % This function performs Bland-Altman analysis on the kinetic energy of
    % two time-resolved, volumetric velocty data sets for all voxels in a provided mask.
    % Returns correlation coefficient, bias, limits of agreement,
    % coefficient of variation, p-value, and vectors to plot the average
    % and difference of KE values.
    
    % Load data
    maskImg = logical(niftiread(maskNiiFile));

    refVelX = niftiread(fullfile(referenceImgNiiDir,'VELX.nii'));
    refVelY = niftiread(fullfile(referenceImgNiiDir,'VELY.nii'));
    refVelZ = niftiread(fullfile(referenceImgNiiDir,'VELZ.nii'));
    testVelX = circshift(niftiread(fullfile(testImgNiiDir,'VELX.nii')),[0,0,0,tShift]);
    testVelY = circshift(niftiread(fullfile(testImgNiiDir,'VELY.nii')),[0,0,0,tShift]);
    testVelZ = circshift(niftiread(fullfile(testImgNiiDir,'VELZ.nii')),[0,0,0,tShift]);
    rho = 1060; % blood density in kg/m^3
    % Compute speed
    speedRef = sqrt(double(refVelX(maskImg)).^2 + double(refVelY(maskImg)).^2 + double(refVelZ(maskImg)).^2);
    speedTest = sqrt(double(testVelX(maskImg)).^2 + double(testVelY(maskImg)).^2 + double(testVelZ(maskImg)).^2);
    % KE density  = 1/2 * rho * v^2
    keRef = 0.5*rho*(speedRef/1000).^2;
    keTest = 0.5*rho*(speedTest/1000).^2;
    
    if (size(refVelX,4)~=size(testVelX,4))
        error('Reference and test images have different number of time frames.');
    end

    diffRT = double(keTest-keRef);
    r = corr(keRef,keTest);
    bias = mean(diffRT);
    loa1 = bias-1.96*std(diffRT);
    loa2 = bias+1.96*std(diffRT);
    CoV = mean(std([keRef,keTest],[],2)./mean([keRef,keTest],2), 'omitnan').*100;
    [~,p] = ttest(keRef-keTest);
    avg = mean([keRef,keTest],2);
end

