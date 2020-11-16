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

function register4DFlow(movingDir,fixedDir,transformFile,outDir)
% This function registers 1 4D flow dataset to another, using
% linear interpolation

    % Load fixed image header
    fInfo = niftiinfo(fullfile(fixedDir,'MAG.nii'));  
        
    % Load moving image
    mInfo = niftiinfo(fullfile(movingDir,'MAG.nii'));
    mMag = niftiread(fullfile(movingDir,'MAG.nii'));
    mVelX = niftiread(fullfile(movingDir,'VELX.nii'));
    mVelY = niftiread(fullfile(movingDir,'VELY.nii'));
    mVelZ = niftiread(fullfile(movingDir,'VELZ.nii'));
    mMagAvg = niftiread(fullfile(movingDir,'AVG_MAG.nii'));
    mVelXAvg = niftiread(fullfile(movingDir,'AVG_VELX.nii'));
    mVelYAvg = niftiread(fullfile(movingDir,'AVG_VELY.nii'));
    mVelZAvg = niftiread(fullfile(movingDir,'AVG_VELZ.nii'));
     
    % Load affine transform
    T1 = loadTransformFile(transformFile);     
    
    % Combine affine transforms
    T0 = (repmat([-1,-1,1,1],4,1).*mInfo.Transform.T);
    T2 = (repmat([-1,-1,1,1],4,1).*fInfo.Transform.T);
    T = T0*pinv(T1)*pinv(T2); % From coord system 1 to coord system 2
    Tinv = pinv(T); % From coord system 2 to coord system 1
    
    % Grab rotation matrix
    R = T(1:3,1:3);
    R = R./repmat(sqrt(sum(R.^2,1)),[3,1]);
    
    % Rearrange velocity directions to line up with fixed image directions
    vX = mVelX*R(1,1)+mVelY*R(2,1)+mVelZ*R(3,1);
    vY = mVelX*R(1,2)+mVelY*R(2,2)+mVelZ*R(3,2);
    vZ = mVelX*R(1,3)+mVelY*R(2,3)+mVelZ*R(3,3);
    vXAvg = mVelXAvg*R(1,1)+mVelYAvg*R(2,1)+mVelZAvg*R(3,1);
    vYAvg = mVelXAvg*R(1,2)+mVelYAvg*R(2,2)+mVelZAvg*R(3,2);
    vZAvg = mVelXAvg*R(1,3)+mVelYAvg*R(2,3)+mVelZAvg*R(3,3);
    
    % Make interpolation grids
    [x1,y1,z1] = ndgrid(0:(mInfo.ImageSize(1)-1),0:(mInfo.ImageSize(2)-1),0:(mInfo.ImageSize(3)-1));
    [x2,y2,z2] = ndgrid(0:(fInfo.ImageSize(1)-1),0:(fInfo.ImageSize(2)-1),0:(fInfo.ImageSize(3)-1));
    xOut = x2*Tinv(1,1) + y2*Tinv(2,1) + z2*Tinv(3,1) + Tinv(4,1);
    yOut = x2*Tinv(1,2) + y2*Tinv(2,2) + z2*Tinv(3,2) + Tinv(4,2);
    zOut = x2*Tinv(1,3) + y2*Tinv(2,3) + z2*Tinv(3,3) + Tinv(4,3);
    
    
    % Interpolate moving images onto fixed image grid
    fMag = int16(zeros(fInfo.ImageSize)); fVelX=fMag; fVelY=fMag; fVelZ=fMag;
    for tt = 1:fInfo.ImageSize(4)
        fMag(:,:,:,tt) = int16(interp3(y1,x1,z1,double(mMag(:,:,:,tt)),yOut,xOut,zOut, 'nearest', 0));
        fVelX(:,:,:,tt) = int16(interp3(y1,x1,z1,double(vX(:,:,:,tt)),yOut,xOut,zOut, 'linear', 0));
        fVelY(:,:,:,tt) = int16(interp3(y1,x1,z1,double(vY(:,:,:,tt)),yOut,xOut,zOut, 'linear', 0));
        fVelZ(:,:,:,tt) = int16(interp3(y1,x1,z1,double(vZ(:,:,:,tt)),yOut,xOut,zOut, 'linear', 0));
    end
    fMagAvg = int16(interp3(y1,x1,z1,double(mMagAvg),yOut,xOut,zOut, 'nearest', 0));
    fVelXAvg = int16(interp3(y1,x1,z1,double(vXAvg),yOut,xOut,zOut, 'linear', 0));
    fVelYAvg = int16(interp3(y1,x1,z1,double(vYAvg),yOut,xOut,zOut, 'linear', 0));
    fVelZAvg = int16(interp3(y1,x1,z1,double(vZAvg),yOut,xOut,zOut, 'linear', 0));
    
    % Write registered images
    if ~exist(outDir,'dir'); mkdir(outDir); end
    niftiwrite(fMag,fullfile(outDir,'MAG.nii'),fInfo);
    niftiwrite(fVelX,fullfile(outDir,'VELX.nii'),fInfo);
    niftiwrite(fVelY,fullfile(outDir,'VELY.nii'),fInfo);
    niftiwrite(fVelZ,fullfile(outDir,'VELZ.nii'),fInfo);
    
    fInfo.ImageSize = fInfo.ImageSize(1:3);
    fInfo.PixelDimensions = fInfo.PixelDimensions(1:3);
    fInfo.TimeUnits = 'None';
        
    niftiwrite(fMagAvg,fullfile(outDir,'AVG_MAG.nii'),fInfo);
    niftiwrite(fVelXAvg,fullfile(outDir,'AVG_VELX.nii'),fInfo);
    niftiwrite(fVelYAvg,fullfile(outDir,'AVG_VELY.nii'),fInfo);
    niftiwrite(fVelZAvg,fullfile(outDir,'AVG_VELZ.nii'),fInfo);
end

% Load an affine transform text file produced in ITK-snap
function T = loadTransformFile(file)
    fid = fopen(file);
    fgetl(fid); fgetl(fid); fgetl(fid);
    line = fgetl(fid);
    fclose(fid);
    line = split(line);
    line = line(2:end);
    line = cellfun(@str2num,line);
    line = reshape(line,[3,4]);
    T = [line;[0,0,0,1]]'; 
    T(1:3,1:3) = T(1:3,1:3)';
end