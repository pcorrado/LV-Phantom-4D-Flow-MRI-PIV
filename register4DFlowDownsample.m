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

function register4DFlowDownsample(movingDir,fixedDir,transformFile,outDir)
% This function registers 1 4D flow dataset to another, where the fixed
% image is at a lower resolution than the moving imageset, using a gaussian
% gridding kernel to grid the values of the moving image to the coordinate
% system of the fixed image.

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
    
    % Compute gaussian kernel and pad moving image
    sigmaX = fInfo.PixelDimensions(1)/mInfo.PixelDimensions(1);
    sigmaY = fInfo.PixelDimensions(2)/mInfo.PixelDimensions(2);
    sigmaZ = fInfo.PixelDimensions(3)/mInfo.PixelDimensions(3);
    [xxx,yyy,zzz] = ndgrid(round(-1.5*sigmaX):round(1.5*sigmaX),round(-1.5*sigmaY):round(1.5*sigmaY),round(-1.5*sigmaZ):round(1.5*sigmaZ));
    gaussKernel = exp( -( (xxx.^2)/(2*sigmaX) + (yyy.^2)/(2*sigmaY) + (zzz.^2)/(2*sigmaZ) ) );
    gaussKernel = gaussKernel./sum(gaussKernel(:));
    padX = max(abs(xxx(:)));
    padY = max(abs(xxx(:)));
    padZ = max(abs(xxx(:)));
    mMag = padarray(mMag,[padX,padY,padZ],0,'both');
    vX = padarray(vX,[padX,padY,padZ],0,'both');
    vY = padarray(vY,[padX,padY,padZ],0,'both');
    vZ = padarray(vZ,[padX,padY,padZ],0,'both');
    
    % Interpolate moving images onto fixed image grid
    fMag = zeros(fInfo.ImageSize); fVelX=fMag; fVelY=fMag; fVelZ=fMag;
    for tt = 1:fInfo.ImageSize(4)
        fprintf('T = %i\n',tt); 
        for zz = 1:fInfo.ImageSize(3)
            for yy = 1:fInfo.ImageSize(2)
                for xx = 1:fInfo.ImageSize(1)
                    % Mu is the center posistion of the gridding kernel in
                    % the moving image
                    muX = xOut(xx,yy,zz);
                    muY = yOut(xx,yy,zz);
                    muZ = zOut(xx,yy,zz);
                    
                    % If mu lands within the image, perform gridding
                    if muX>=1 && muX<=mInfo.ImageSize(1) && muY>=1 && muY<=mInfo.ImageSize(2) && muZ>=1 && muZ<=mInfo.ImageSize(3)
                        xxxx = round((round(-1.5*sigmaX):round(1.5*sigmaX)) + muX + padX);
                        yyyy = round((round(-1.5*sigmaY):round(1.5*sigmaY)) + muY + padY);
                        zzzz = round((round(-1.5*sigmaZ):round(1.5*sigmaZ)) + muZ + padZ);
                        
                        fMag(xx,yy,zz,tt) = sum(sum(sum(double(mMag(xxxx,yyyy,zzzz,tt)).*gaussKernel,1),2),3);
                        fVelX(xx,yy,zz,tt) = sum(sum(sum(double(vX(xxxx,yyyy,zzzz,tt)).*gaussKernel,1),2),3);
                        fVelY(xx,yy,zz,tt) = sum(sum(sum(double(vY(xxxx,yyyy,zzzz,tt)).*gaussKernel,1),2),3);
                        fVelZ(xx,yy,zz,tt) = sum(sum(sum(double(vZ(xxxx,yyyy,zzzz,tt)).*gaussKernel,1),2),3);
                    end
                end
            end
        end
    end
    % Use linear interpolation for the average images (these are not used
    % in manuscript).
    fMagAvg = int16(interp3(y1,x1,z1,double(mMagAvg),yOut,xOut,zOut, 'nearest', 0));
    fVelXAvg = int16(interp3(y1,x1,z1,double(vXAvg),yOut,xOut,zOut, 'linear', 0));
    fVelYAvg = int16(interp3(y1,x1,z1,double(vYAvg),yOut,xOut,zOut, 'linear', 0));
    fVelZAvg = int16(interp3(y1,x1,z1,double(vZAvg),yOut,xOut,zOut, 'linear', 0));
    

    
    % Write registered images
    if ~exist(outDir,'dir'); mkdir(outDir); end
    niftiwrite(int16(round(fMag)),fullfile(outDir,'MAG.nii'),fInfo);
    niftiwrite(int16(round(fVelX)),fullfile(outDir,'VELX.nii'),fInfo);
    niftiwrite(int16(round(fVelY)),fullfile(outDir,'VELY.nii'),fInfo);
    niftiwrite(int16(round(fVelZ)),fullfile(outDir,'VELZ.nii'),fInfo);
    
    fInfo.ImageSize = fInfo.ImageSize(1:3);
    fInfo.PixelDimensions = fInfo.PixelDimensions(1:3);
    fInfo.TimeUnits = 'None';
        
    niftiwrite(int16(round(fMagAvg)),fullfile(outDir,'AVG_MAG.nii'),fInfo);
    niftiwrite(int16(round(fVelXAvg)),fullfile(outDir,'AVG_VELX.nii'),fInfo);
    niftiwrite(int16(round(fVelYAvg)),fullfile(outDir,'AVG_VELY.nii'),fInfo);
    niftiwrite(int16(round(fVelZAvg)),fullfile(outDir,'AVG_VELZ.nii'),fInfo);
end

% Load affine transform matrix from ITK-snap
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