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

function KE = computeKineticEnergy(maskFile,imgDir)
% This is a function to compute the average kinetic energy in microJoules 
% of a time-resolved, 3D velocity map in a given mask.

    % Load images
    mask = logical(niftiread(maskFile));
    rho = 1060; % blood density in kg/m^3
    
    info = niftiinfo(fullfile(imgDir,'VELX.nii'));
    velX = niftiread(fullfile(imgDir,'VELX.nii'));
    velY = niftiread(fullfile(imgDir,'VELY.nii'));
    velZ = niftiread(fullfile(imgDir,'VELZ.nii'));
    vol = prod(info.PixelDimensions(1:3)/1000); % pixel volume in m^3

    % Loop through time frames
    for t=1:size(velZ,4)
        m = mask(:,:,:,t);

        vX = double(velX(:,:,:,t))./1000; % in m/s
        vY = double(velY(:,:,:,t))./1000; % in m/s
        vZ = double(velZ(:,:,:,t))./1000; % in m/s
        % Compute speed
        vMag = sqrt(vX.^2 + vY.^2 + vZ.^2); % in m/s
        
        % Compute kinetic energy
        KE(t) = (1/2)*rho*vol*sum(vMag(m).^2).*1e6; %#ok<AGROW> % in uJ
    end
    KE = mean(KE);
end