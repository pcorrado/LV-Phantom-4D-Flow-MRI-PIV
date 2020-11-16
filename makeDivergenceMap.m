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

function makeDivergenceMap(srcDir)
% This function computes divergence voxel-by-voxel on MR or PIV images.
% div(V) = dVx/dx + dVy/dy + dVz/dz
    
   curDir = pwd; 
    if nargin<1 || (~ischar(srcDir) && ~isstring(srcDir)) || ~exist(srcDir, 'dir')
        srcDir = curDir;    
        fprintf('Setting source directory to %s\n', srcDir);
    end
    srcDir = char(srcDir);
    
    info = niftiinfo(fullfile(srcDir,'VELX.nii'));
    % Load and negate data because my MRI data is stored so that dark 
    % voxels denote positive velocities
    velX = -niftiread(fullfile(srcDir,'VELX.nii'));
    velY = -niftiread(fullfile(srcDir,'VELY.nii'));
    velZ = -niftiread(fullfile(srcDir,'VELZ.nii'));
    
    [X, Y, Z] = ndgrid((1:info.ImageSize(1)).*info.PixelDimensions(1),...
                         (1:info.ImageSize(2)).*info.PixelDimensions(2),...
                         (1:info.ImageSize(3)).*info.PixelDimensions(3));
      
    % Compute discrete derivatives
    dx = circshift(X,[1,0,0,0])-X;
    dy = circshift(Y,[0,1,0,0])-Y;
    dz = circshift(Z,[0,0,1,0])-Z;
    dVx = circshift(velX,[1,0,0,0])-velX;
    dVy = circshift(velY,[0,1,0,0])-velY;
    dVz = circshift(velZ,[0,0,1,0])-velZ;
    
    % Compute divergence
    div = int16(double(dVx)./dx + double(dVy)./dy + double(dVz)./dz);
    
    
    info.Datatype = class(div);
    niftiwrite(div,fullfile(srcDir,'DIVERGENCE.nii'),info,'Compressed',false);
    
end