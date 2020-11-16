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

% This script registers the Day 1 MRI data to Day 2 MRI data in order to
% perform repeatablity comparisons
      
%update this to match your file structure
MRI_Dirs = ["E:\MRI_Data\Day 1\PC VIPR",...
            "E:\MRI_Data\Day 1\PC SOS";...
            "E:\MRI_Data\Day 2\PC VIPR",...
            "E:\MRI_Data\Day 2\PC SOS"];
transform = "E:\MRI_Data\MRIDay1_to_MRIDay2_Transform.txt";
        
trajectories = ["PC VIPR","PC SOS"];
times = ["30_Minutes","20_Minutes","10_Minutes","5_Minutes","2.5_Minutes"];

% Loop through trajectories
for d = 1:size(MRI_Dirs,1)
    % Loop through scan times
    for t = 1:numel(times)
        fprintf('Registering %s.\n', fullfile(MRI_Dirs(1,d),times(t)));
        % Change the output directory as appropriate
        register4DFlow(fullfile(MRI_Dirs(1,d),times(t)),...
                       fullfile(MRI_Dirs(2,d),times(t)),...
                       transform,...
                       sprintf('E:\\MRI_Data\\Day1RegisteredToDay2\\%s\\%s',trajectories(d),times(t)));
    end
end