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

% This script was used to register the PIV data to the MRI images 

averageAll();
PIV_DIRS = [...
            "E:\PIV_Data\Unregistered\Heartbeat 1";...
            "E:\PIV_Data\Unregistered\Heartbeat 2";...
            "E:\PIV_Data\Unregistered\Heartbeat 3";...
            "E:\PIV_Data\Unregistered\Heartbeat 4";...
            "E:\PIV_Data\Unregistered\Average Beat";...
            ];
        
MRI_DIRS = ["E:\MRI_Data\Day 1\PC VIPR\30_Minutes",...
            "E:\MRI_Data\Day 1\PC SOS\30_Minutes";...
            "E:\MRI_Data\Day 2\PC VIPR\30_Minutes",...
            "E:\MRI_Data\Day 2\PC SOS\30_Minutes"];

TRANSFORMS = ["E:\PIV_Data\PIV2MRI_Transform_Day1.txt","E:\PIV_Data\PIV2MRI_Transform_Day2.txt"];
BASEOUTPUT = "E:\PIV_Data";

days = ["Registered to Day 1"; "Registered to Day 2"];
trajectories = ["PC VIPR", "PC SOS"];

beats = [...
    "Heartbeat 1";...
    "Heartbeat 2";...
    "Heartbeat 3";...
    "Heartbeat 4";...
    "Average Beat";...
    ];

% Loop through days and trajectories and register PIV set to each MRI
% geometry, using gaussian gridding kernel.
for d = 1:numel(days)
    for t = 1:numel(trajectories)
        for hb = 1:numel(PIV_DIRS)
            register4DFlowDownsample(PIV_DIRS(hb),MRI_DIRS(d,t),TRANSFORMS(d),sprintf('%s\\%s\\%s\\%s',BASEOUTPUT, days(d),trajectories(t),beats(hb)));
        end
    end
end

function averageAll() % Combine heartbeats into average heartbeat
    averageNii('MAG.nii');
    averageNii('VELX.nii');
    averageNii('VELY.nii');
    averageNii('VELZ.nii');
    averageNii('AVG_MAG.nii');
    averageNii('AVG_VELX.nii');
    averageNii('AVG_VELY.nii');
    averageNii('AVG_VELZ.nii');
end
function averageNii(fileName) % Combine file from all heartbeats into average heartbeat
    PIV_Dirs = [...
            "E:\PIV_Raw\YesShift\Unregistered\Heartbeat 1";...
            "E:\PIV_Raw\YesShift\Unregistered\Heartbeat 2";...
            "E:\PIV_Raw\YesShift\Unregistered\Heartbeat 3";...
            "E:\PIV_Raw\YesShift\Unregistered\Heartbeat 4"];
        
    info = niftiinfo(fullfile(PIV_Dirs(1),fileName));
    img = niftiread(fullfile(PIV_Dirs(1),fileName));
    img = img+niftiread(fullfile(PIV_Dirs(2),fileName));
    img = img+niftiread(fullfile(PIV_Dirs(3),fileName));
    img = img+niftiread(fullfile(PIV_Dirs(4),fileName));
    img = int16(img./4);
    niftiwrite(img,fullfile("E:\PIV_Raw\YesShift\Unregistered\Average Beat",fileName),info);
end