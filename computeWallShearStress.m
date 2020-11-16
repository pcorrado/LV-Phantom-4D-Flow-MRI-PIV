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

function wss = computeWallShearStress(maskFile,imgDir,smoothMaskMeshFile)
% Compute average volumetric wall shear stress according to the derivation 
% provided in the appendix of Stalder AF, et al. Quantitative 2D and 3D phase contrast MRI: 
% Optimized analysis of blood flow and vessel wall parameters. Magn Reson 
% Med. 2008;60:1218–1231. doi: 10.1002/mrm.21778
% However, since we are computing wall shear stress along a whole surface
% rather than just at a 2D plane, we do not use Green's Theorem and rather
% just compute the wall shear stress vector based on discrete derivatives
% for each voxel touching the surface.
%
% WSS = [2*n1*dVx/dx + n2*(dVy/dx dVx/dy) + n3*(dVz/dx + dVx/dz);
%        n1*(dVy/dx + dVx/dy) + 2*n2*dVy/dy + n3*(dVz/dy + dVy/dz); 
%        n1*(dVz/dx + dVx/dz) + n2*(dVy/dz dVz/dy) + 2*n3*dVz/dz]

    % Load data
    info = niftiinfo(fullfile(imgDir,'VELX.nii'));
    velX = niftiread(fullfile(imgDir,'VELX.nii'));
    velY = niftiread(fullfile(imgDir,'VELY.nii'));
    velZ = niftiread(fullfile(imgDir,'VELZ.nii'));

    nu = 0.00824; % viscosity in Pa*sec
    
    % Get location and surface normal for each voxel touching the surface
    % mesh.
    [n,ind] = getSmoothMask(maskFile,smoothMaskMeshFile, info.ImageSize, info.PixelDimensions);
      
    wss = zeros(size(velZ,4),1);
    
    % Loop through time frames
    for t=1:size(velZ,4)
        vX = double(velX(:,:,:,t))./1000; % in m/s
        vY = double(velY(:,:,:,t))./1000; % in m/s
        vZ = double(velZ(:,:,:,t))./1000; % in m/s
        
        % take discrete derivatives in units of m/s/m
        [dVxdy,dVxdx,dVxdz] = gradient(vX./(info.PixelDimensions(1)/1000));
        [dVydy,dVydx,dVydz] = gradient(vY./(info.PixelDimensions(2)/1000));
        [dVzdy,dVzdx,dVzdz] = gradient(vZ./(info.PixelDimensions(3)/1000));
        
        % Sample gradients at surface voxels and put into cell matrix
        dVdx = {dVxdx(ind),dVxdy(ind),dVxdz(ind); dVydx(ind),dVydy(ind),dVydz(ind); dVzdx(ind),dVzdy(ind),dVzdz(ind)};
        wss = zeros(size(n,1),3);
        %Loop over velocity directions
        for ii=1:3
           wss(:,ii) = 2.*nu.*n(:,ii).*dVdx{ii,ii};
           others = setdiff(1:3,ii);
           for iii=1:2
               jj = others(iii);
               wss(:,ii) = wss(:,ii) + nu.*n(:,jj).*(dVdx{ii,jj}+dVdx{jj,ii});
           end
        end
        
        % Take magnitude of the wss vector
        wss = sqrt(wss(:,1).^2+wss(:,2).^2+wss(:,3).^2);
        % Take the average at each time point
        wss(t) = mean(wss);
    end
    wss = mean(wss)*1000; %Take mean and convert from Pa to mPa
end

% Function to creat smooth surface mesh with top cut off for LV WSS
% computation
function [n,ind] = smoothMask(maskFile,smoothMaskFile, imgSize, pixDims)
    mask = squeeze(niftiread(maskFile));
    [x,y,z] = ndgrid((1:imgSize(1))*pixDims(1),(1:imgSize(2))*pixDims(2),(1:imgSize(3))*pixDims(3));
    lvInd = mask(:,:,:,1)>0;

    mesh = alphaShape(x(lvInd),y(lvInd),z(lvInd));
    [tri, xyz] = boundaryFacets(mesh);
    FV.faces= tri;
    FV.vertices=xyz;
    FV2 = smoothpatch(FV);
    tri=FV2.faces;
    xyz = FV2.vertices;
   
    [allZ,ind] = sort((xyz(tri(:,1),3)+xyz(tri(:,2),3)+xyz(tri(:,3),3))/3);
    last = find((allZ-min(allZ))/(max(allZ)-min(allZ))>0.97,1);
    tri = tri(ind(1:last),:);
    
%     trisurf(tri,xyz(:,1),xyz(:,2),xyz(:,3),'FaceColor','cyan','FaceAlpha',0.3);

    X = mean([xyz(tri(:,1),1),xyz(tri(:,2),1),xyz(tri(:,3),1)],2);
    Y = mean([xyz(tri(:,1),2),xyz(tri(:,2),2),xyz(tri(:,3),2)],2);
    Z = mean([xyz(tri(:,1),3),xyz(tri(:,2),3),xyz(tri(:,3),3)],2);
    XX=round(X./pixDims(1));
    YY=round(Y./pixDims(2));
    ZZ=round(Z./pixDims(3));
    d1 = [xyz(tri(:,2),1)-xyz(tri(:,1),1),...
          xyz(tri(:,2),2)-xyz(tri(:,1),2),...
          xyz(tri(:,2),3)-xyz(tri(:,1),3),];
    d2 = [xyz(tri(:,3),1)-xyz(tri(:,1),1),...
          xyz(tri(:,3),2)-xyz(tri(:,1),2),...
          xyz(tri(:,3),3)-xyz(tri(:,1),3),];
    n = cross(d1,d2);
    n = n./repmat(sqrt(n(:,1).^2 + n(:,2).^2 + n(:,3).^2),1,3);
    ind = sub2ind(imgSize(1:3),XX,YY,ZZ);
    save(smoothMaskFile,'n','ind');
end

% Function to load surface mesh
function [n,ind] = loadSmoothMask(smoothMaskFile)
    dat=load(smoothMaskFile);
    n=dat.n;
    ind=dat.ind;
end

% If mesh already exists, load it, if not, create one.
function [n,ind] = getSmoothMask(maskFile,smoothMaskFile, imgSize, pixDims)
    if exist(smoothMaskFile,'file')
        [n,ind] = loadSmoothMask(smoothMaskFile);
    else
        [n,ind] = smoothMask(maskFile,smoothMaskFile, imgSize, pixDims);
    end
end