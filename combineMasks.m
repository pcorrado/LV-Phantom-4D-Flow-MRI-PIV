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

% This is a script for combining masks segmented both on PIV data and MRI
% data, and keeping as a left ventricle mask the intersection of the 2
% masks. 
clear all; close all;  %#ok<CLALL>

PIV_Dirs = [...
            "E:\PIV Data\Registered to Day 1\PC VIPR\Average Beat";...
            "E:\PIV Data\Registered to Day 1\PC SOS\Average Beat";...
            "E:\PIV Data\Registered to Day 2\PC VIPR\Average Beat";...
            "E:\PIV Data\Registered to Day 2\PC SOS\Average Beat";...
            ];
        
MRI_masks = ["E:\MRI_Data\Day 1\PC VIPR\Mask_VIPR_Day1.nii.gz";...
             "E:\MRI_Data\Day 1\PC SOS\Mask_SOS_Day1.nii.gz";...
             "E:\MRI_Data\Day 2\PC VIPR\Mask_VIPR_Day2.nii.gz";...
             "E:\MRI_Data\Day 2\PC SOS\Mask_SOS_Day2.nii.gz"];

for ii = 1:numel(PIV_Dirs)
    magFiles = dir(fullfile(PIV_Dirs(ii),'MAG.nii'));
    for iii=1:numel(magFiles)
        mag = niftiread(fullfile(PIV_Dirs(ii),'MAG.nii'));
        magInfo = niftiinfo(fullfile(PIV_Dirs(ii),'MAG.nii'));
        mask = niftiread(MRI_masks(ii));
        mask = mag.*int16(mask);
        niftiwrite(mask,fullfile(PIV_Dirs(ii),'..\MASK.nii'),magInfo);
    end
end