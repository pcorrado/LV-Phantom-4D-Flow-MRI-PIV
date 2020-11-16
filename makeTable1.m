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

% This script reproduces Table 1 from the manuscript, which
% tabulates divergence, speed, kinetic energy, and wall shear stress data
% for the averge PIV dataset plus all MRI datasets.

clear all; close all; %#ok<CLALL>

% UPDATE THIS BASED ON YOUR FILE STRUCTURE
BASE_PIV_DIR = "E:\PIV_Data";
BASE_MRI_DIR = "E:\MRI_Data";

piv_days = ["Registered to Day 1"; "Registered to Day 2"];
trajectories = ["PC VIPR", "PC SOS"];

scanTime = ["30","20","10","5","2.5"];
mri_days = ["Day 1", "Day 2"];

R = [14	21	41	83	165	5	8	17	33	66];

T = table('Size',[20,12],'VariableTypes',{'string','string','string','string','string','string','string','string','string','string','string','string'},...
'VariableNames',{'Acquisition','PIV','PCVIPR_30','PCVIPR_20','PCVIPR_10','PCVIPR_5','PCVIPR_2p5','PCSOS_30','PCSOS_20','PCSOS_10','PCSOS_5','PCSOS_2p5'});

T(:,1) = {...
    'R';...
    'Scan time (mins)';...
    'RMS Div. ((cm/s)/cm)';...
    'Speed Corr. Coeff.';...
    'Speed Bias (cm/s)';...
    'Speed LoA (cm/s)';...
    'Speed CoV (%)';...
    'Speed P-Value';...
    'Mean KE (µJ)';...
    'KE Corr. Coeff.';...
    'KE Bias (nJ/mL)';...
    'KE LoA (nJ/mL)';...
    'KE CoV (%)';...
    'KE P-Value';...
    'Mean WSS (mPa)';...
    'WSS Corr. Coeff.';...
    'WSS Bias (mPa)';...
    'WSS LoA (mPa)';...
    'WSS CoV (%)';...
    'WSS P-Value'};

% Loop through all image sets and fill in table
for col=2:(2+numel(trajectories)*numel(scanTime))
    fprintf('Analyzing column %i of %i.\n',col-1,1+numel(trajectories)*numel(scanTime));
    
    if col==2 % if PIV dataset
        traj='PIV';
        pivDir = string(sprintf('%s\\%s\\%s\\Average Beat',BASE_PIV_DIR,piv_days(1),trajectories(1)));
        maskFile = string(fullfile(pivDir,'..\MASK.nii'));
        wallMesh = "C:\\Users\\pcorrado\\Desktop\\smoothMask_PC VIPR_Day1.mat";

        % Since other metrics are computed relative to PIV, PIV column only
        % has results for mean divergence, KE, and WSS.
        T(1,col) = {'--'};
        T(2,col) = {'--'};   
        T(3,col) = {sprintf('%2.2f',computeDivergence(maskFile,pivDir))};
        
        T(4,col) = {'--'};
        T(5,col) = {'--'};
        T(6,col) = {'--'};
        T(7,col) = {'--'};
        T(8,col) = {'--'};
        
        T(9,col) = {sprintf('%2.0f',computeKineticEnergy(maskFile,pivDir))};
        T(10,col) = {'--'};
        T(11,col) = {'--'};
        T(12,col) = {'--'};
        T(13,col) = {'--'};
        T(14,col) = {'--'};
        
        T(15,col) = {sprintf('%2.2f',computeWallShearStress(maskFile,pivDir,wallMesh))};
        T(16,col) = {'--'};
        T(17,col) = {'--'};
        T(18,col) = {'--'};
        T(19,col) = {'--'};
        T(20,col) = {'--'};
    else
        trajNum = ceil((col-2)/numel(scanTime));
        traj = trajectories(trajNum);
        % MRI day 1
        pivDir(1) = string(sprintf('%s\\%s\\%s\\Average Beat',BASE_PIV_DIR,piv_days(1),trajectories(trajNum)));
        maskFile(1) = string(fullfile(pivDir(1),'..\MASK.nii'));
        wallMesh(1) = string(sprintf('C:\\Users\\pcorrado\\Desktop\\smoothMask_%s_Day1.mat',traj));
        mriDir(1) = string(sprintf('%s\\%s\\%s\\%s_Minutes\\',BASE_MRI_DIR,mri_days(1),traj,scanTime(mod(col-3,numel(scanTime))+1)));
        % MRI day 2
        pivDir(2) = string(sprintf('%s\\%s\\%s\\Average Beat',BASE_PIV_DIR,piv_days(2),trajectories(trajNum)));
        maskFile(2) = string(fullfile(pivDir(2),'..\MASK.nii'));
        wallMesh(2) = string(sprintf('C:\\Users\\pcorrado\\Desktop\\smoothMask_%s_Day2.mat',traj));
        mriDir(2) = string(sprintf('%s\\%s\\%s\\%s_Minutes\\',BASE_MRI_DIR,mri_days(2),traj,scanTime(mod(col-3,numel(scanTime))+1)));
        
        % Acceleration factor, scan time, and diverence
        T(1,col) = {R(col-2)};
        T(2,col) = {scanTime(mod(col-3,numel(scanTime))+1)};
        T(3,col) = {getMeanStandardDeviation(@computeDivergence,maskFile(1),mriDir(1),maskFile(2),mriDir(2))};

        % Speed metrics
        [r(1),bias(1),loa1(1),loa2(1),CoV(1),p(1),~,~] = computeSpeedBlandAltman(maskFile(1),pivDir(1),mriDir(1),0);
        [r(2),bias(2),loa1(2),loa2(2),CoV(2),p(2),~,~] = computeSpeedBlandAltman(maskFile(2),pivDir(2),mriDir(2),0);
        T(4,col) = {sprintf('%2.2f (%2.2f)',mean(r),std(r))};
        T(5,col) = {sprintf('%1.1f (%1.1f)',mean(bias)/10,std(bias)/10)};
        T(6,col) = {sprintf('[%1.1f, %1.1f]',mean(loa1)/10,mean(loa2)/10)};
        T(7,col) = {sprintf('%2.0f (%1.0f)',mean(CoV),std(CoV))};
        T(8,col) = {sprintf('%3.3f',mean(p))};
        
        % KE Metrics
        [r(1),bias(1),loa1(1),loa2(1),CoV(1),p(1),~,~] = computeKEBlandAltman(maskFile(1),pivDir(1),mriDir(1),0);
        [r(2),bias(2),loa1(2),loa2(2),CoV(2),p(2),~,~] = computeKEBlandAltman(maskFile(2),pivDir(2),mriDir(2),0);
        T(9,col) = {getMeanStandardDeviation(@computeKineticEnergy,maskFile(1),mriDir(1),maskFile(2),mriDir(2))};
        T(10,col) = {sprintf('%2.2f (%2.2f)',mean(r),std(r))};
        T(11,col) = {sprintf('%1.1f (%1.1f)',mean(bias),std(bias))};
        T(12,col) = {sprintf('[%1.1f, %1.1f]',mean(loa1),mean(loa2))};
        T(13,col) = {sprintf('%2.0f (%1.0f)',mean(CoV),std(CoV))};
        T(14,col) = {sprintf('%3.3f',mean(p))};
        
        % WSS Metrics
        T(15,col) = {getMeanStandardDeviation(@computeWallShearStress,maskFile(1),mriDir(1),wallMesh(1),maskFile(2),mriDir(2),wallMesh(2))};
        [r(1),bias(1),loa1(1),loa2(1),CoV(1),p(1)] = computeWSSBlandAltman(maskFile(1),pivDir(1),mriDir(1),0,wallMesh(1));
        [r(2),bias(2),loa1(2),loa2(2),CoV(2),p(2)] = computeWSSBlandAltman(maskFile(2),pivDir(2),mriDir(2),0,wallMesh(2));
        T(16,col) = {sprintf('%2.2f (%2.2f)',mean(r),std(r))};
        T(17,col) = {sprintf('%1.1f (%1.1f)',mean(bias),std(bias))};
        T(18,col) = {sprintf('[%1.1f, %1.1f]',mean(loa1),mean(loa2))};
        T(19,col) = {sprintf('%2.0f (%1.0f)',mean(CoV),std(CoV))};
        T(20,col) = {sprintf('%3.3f',mean(p))};
    end
    disp(T);
end

% Format data for table so you compute metric for Day 1 and Day 2, then
% report the mean and the standard deviation
function meanAndStd = getMeanStandardDeviation(varargin)
    fnctn = varargin{1};
    numArgs = (nargin-1)/2;
    val1 = fnctn(varargin{2:(numArgs+1)});
    val2 = fnctn(varargin{(numArgs+2):end});
    meanAndStd = sprintf('%2.2f (%2.2f)',mean([val1,val2]),std([val1,val2]));
end

% Get mean divergence within mask from divergence map
function divVal = computeDivergence(maskFile,imgDir)
    divFile = fullfile(imgDir,'DIVERGENCE.nii');
    if ~exist(divFile,'file')
        makeDivergenceMap(imgDir);
    end
    divSquared = double(niftiread(divFile)).^2;
    maskImg = logical(niftiread(maskFile));
    divVal = sqrt(mean(divSquared(maskImg~=0)))/10;
end