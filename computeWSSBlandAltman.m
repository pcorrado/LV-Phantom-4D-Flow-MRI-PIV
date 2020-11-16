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

function [r,bias,loa1,loa2,CoV,p] = computeWSSBlandAltman(maskFile, referenceImgNiiDir, testImgNiiDir, tshift, smoothMaskMeshFile)
% This function performs Bland-Altman analysis on the wall shear stress of 
% two time-resolved, volumetric velocty data sets for all voxels in a provided mask.
% Returns correlation coefficient, bias, limits of agreement,
% coefficient of variation, p-value, and vectors to plot the average
% and difference of WSS values.
% Wall shear stress  computation uses the derivation provided in the 
% appendix of Stalder AF, et al. Quantitative 2D and 3D phase contrast MRI: 
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
    info = niftiinfo(fullfile(referenceImgNiiDir,'VELX.nii'));
    refVelX = niftiread(fullfile(referenceImgNiiDir,'VELX.nii'));
    refVelY = niftiread(fullfile(referenceImgNiiDir,'VELY.nii'));
    refVelZ = niftiread(fullfile(referenceImgNiiDir,'VELZ.nii'));
    testVelX = circshift(niftiread(fullfile(testImgNiiDir,'VELX.nii')),[0,0,0,tshift]);
    testVelY = circshift(niftiread(fullfile(testImgNiiDir,'VELY.nii')),[0,0,0,tshift]);
    testVelZ = circshift(niftiread(fullfile(testImgNiiDir,'VELZ.nii')),[0,0,0,tshift]);

    % Get location and normal of each triange on surface mesh of LV
    % segmentation.
    [n,ind] = getSmoothMask(maskFile,smoothMaskMeshFile, info.ImageSize, info.PixelDimensions);
      
    wssRef=[];
    wssTest = [];
    for t=1:size(refVelZ,4)
        vXR = double(refVelX(:,:,:,t))./1000; % in m/s
        vYR = double(refVelY(:,:,:,t))./1000; % in m/s
        vZR = double(refVelZ(:,:,:,t))./1000; % in m/s
        vXT = double(testVelX(:,:,:,t))./1000; % in m/s
        vYT = double(testVelY(:,:,:,t))./1000; % in m/s
        vZT = double(testVelZ(:,:,:,t))./1000; % in m/s
        
        wssRef = [wssRef;getWSS(vXR,vYR,vZR,n,ind,info)]; %#ok<AGROW>
        wssTest = [wssTest;getWSS(vXT,vYT,vZT,n,ind,info)]; %#ok<AGROW>
    end
    
    % Compare reference and test datasets 
    diffRT = double(wssTest-wssRef);
    r = corr(wssRef,wssTest);
    bias = mean(diffRT);
    loa1 = bias-1.96*std(diffRT);
    loa2 = bias+1.96*std(diffRT);
    CoV = mean(std([wssRef,wssTest],[],2)./mean([wssRef,wssTest],2), 'omitnan').*100;
    [~,p] = ttest(wssRef-wssTest);
end

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

function [n,ind] = loadSmoothMask(smoothMaskFile)
    dat=load(smoothMaskFile);
    n=dat.n;
    ind=dat.ind;
end

function [n,ind] = getSmoothMask(maskFile,smoothMaskFile, imgSize, pixDims)
    if exist(smoothMaskFile,'file')
        [n,ind] = loadSmoothMask(smoothMaskFile);
    else
        [n,ind] = smoothMask(maskFile,smoothMaskFile, imgSize, pixDims);
    end
end

% For comments on how this code works, refer to the function
% computeWallShearStress.m
function wss = getWSS(vX,vY,vZ,n,ind,info)
    nu = 0.00824; % viscosity in Pa*sec

    [dVxdy,dVxdx,dVxdz] = gradient(vX./(info.PixelDimensions(1)/1000));
    [dVydy,dVydx,dVydz] = gradient(vY./(info.PixelDimensions(2)/1000));
    [dVzdy,dVzdx,dVzdz] = gradient(vZ./(info.PixelDimensions(3)/1000));
    
    dVdx = {dVxdx(ind),dVxdy(ind),dVxdz(ind); dVydx(ind),dVydy(ind),dVydz(ind); dVzdx(ind),dVzdy(ind),dVzdz(ind)};
    wss = zeros(size(n,1),3);
    for ii=1:3
        wss(:,ii) = 2.*nu.*n(:,ii).*dVdx{ii,ii};
        others = setdiff(1:3,ii);
        for iii=1:2
            jj = others(iii);
            wss(:,ii) = wss(:,ii) + nu.*n(:,jj).*(dVdx{ii,jj}+dVdx{jj,ii});
        end
    end
    
    wss = sqrt(wss(:,1).^2+wss(:,2).^2+wss(:,3).^2)*1000; % in mPa
end