clear all; close all; %#ok<CLALL>

imgDirs = ["E:\PIV_Raw\YesShift\Registered to Day 1\PC VIPR\Average Beat";...
           "E:\MRI_Data\Day 1\PC VIPR\20_Minutes";...
           "E:\MRI_Data\Day 1\PC SOS\20_Minutes";...
           "E:\MRI_Data\Day 1\PC VIPR\5_Minutes";...
           "E:\MRI_Data\Day 1\PC SOS\5_Minutes"];
       
maskFiles = ["E:\PIV_Raw\YesShift\Registered to Day 1\PC VIPR\MASK.nii";...
             "E:\PIV_Raw\YesShift\Registered to Day 1\PC VIPR\MASK.nii";...
             "E:\PIV_Raw\YesShift\Registered to Day 1\PC SOS\MASK.nii";...
             "E:\PIV_Raw\YesShift\Registered to Day 1\PC VIPR\MASK.nii";...
             "E:\PIV_Raw\YesShift\Registered to Day 1\PC SOS\MASK.nii"];

t = 13;
slice = 26;

imgSize = [126,42];
bigSize = imgSize.*[1,5]+[2,8];
bigImg = ones(bigSize);
 
for col=1:numel(imgDirs)
    fprintf('Analyzing column %i of %i.\n',col,numel(imgDirs));

    mask = niftiread(maskFiles(col));
    pivMask = niftiread(strrep(maskFiles(col),'MASK.nii','Average Beat\MAG.nii'));
%     pivMask = imdilate(pivMask,strel('sphere',1));
    eitherMask = mask + pivMask;
    eitherMask(eitherMask>0) = 1;
    velX = niftiread(fullfile(imgDirs(col),'VELX.nii')).*eitherMask;
    velY = niftiread(fullfile(imgDirs(col),'VELY.nii')).*eitherMask;
    velZ = niftiread(fullfile(imgDirs(col),'VELZ.nii')).*eitherMask;

    velX = fliplr(rot90(squeeze(double(velX(slice,8:49,9:50,t))),1));
    velX = double(velX+150)/300;
    velY = fliplr(rot90(squeeze(double(velY(slice,8:49,9:50,t))),1));
    velY = double(velY+150)/300;
    velZ = fliplr(rot90(squeeze(double(velZ(slice,8:49,9:50,t))),1));
    velZ = double(velZ+150)/300;
    
    mask = fliplr(rot90(squeeze(double(mask(slice,8:49,9:50,t))),1));
    eitherMask = fliplr(rot90(squeeze(double(eitherMask(slice,8:49,9:50,t))),1));
    
    BD = bwboundaries(eitherMask);
    pivContours{col} = BD{1}; %#ok<SAGROW>
    BD = bwboundaries(mask);
    contours{col} = BD{1}; %#ok<SAGROW>
    bigImg(:,((col-1)*(imgSize(2)+2)+1):(col*imgSize(2)+2*(col-1))) = [velX;ones(1,imgSize(2));velY;ones(1,imgSize(2));velZ];
    
end
imshow(bigImg,[0,1]);
hold on;
for col = 1:numel(pivContours)
    plot(pivContours{col}(:,2)+(col-1)*(imgSize(2)+2), pivContours{col}(:,1), '-y', 'LineWidth', 0.75);
    plot(pivContours{col}(:,2)+(col-1)*(imgSize(2)+2), pivContours{col}(:,1)+(imgSize(1)-2)/3+1, '-y', 'LineWidth', 0.75);
    plot(pivContours{col}(:,2)+(col-1)*(imgSize(2)+2), pivContours{col}(:,1)+((imgSize(1)-2)*2)/3+3, '-y', 'LineWidth', 0.75);
    
    plot(contours{col}(:,2)+(col-1)*(imgSize(2)+2), contours{col}(:,1), '--r', 'LineWidth', 0.75);
    plot(contours{col}(:,2)+(col-1)*(imgSize(2)+2), contours{col}(:,1)+(imgSize(1)-2)/3+1, '--r', 'LineWidth', 0.75);
    plot(contours{col}(:,2)+(col-1)*(imgSize(2)+2), contours{col}(:,1)+((imgSize(1)-2)*2)/3+3, '--r', 'LineWidth', 0.75);
end
hold off;
