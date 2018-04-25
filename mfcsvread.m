function outputStuff = mfcsvread(fileName)
%mfcsvread reads a CSV file containing both text & numeric data.  MATLAB's
%csvread function will work with all numeric, or all text data, but not
%both.  It's common to have a file with a single line of comma separated
%text headers followed by many rows of numeric data.  xlsread is limited in
%the number of rows & colums (actually, Excel is the limitation) it can
%read.
%
% The CSV file should look like:
% comma, separated, text, ...
% 1,2,3,4,5,...
% 6,7,8,9,10,...
% etc...
%
% The output is a structure with the column headers as fields in a
% structure each with a vector of data.


% % Copyright (c) 2010, The MathWorks, Inc.
% % All rights reserved.
% % 
% % Redistribution and use in source and binary forms, with or without 
% % modification, are permitted provided that the following conditions are 
% % met:
% % 
% %     * Redistributions of source code must retain the above copyright 
% %       notice, this list of conditions and the following disclaimer.
% %     * Redistributions in binary form must reproduce the above copyright 
% %       notice, this list of conditions and the following disclaimer in 
% %       the documentation and/or other materials provided with the distribution
% %       
% % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% % ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% % LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% % SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% % INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% % CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% % ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% % POSSIBILITY OF SUCH DAMAGE.

% Open the file
fid=fopen(fileName);

% Start reading the data
tline = fgetl(fid); % Read in the first line only (text headers)
tline = tline(tline~=' '); %Get rid of any spaces
commaLocs=findstr(',',tline); % find the commas
fieldNames=cell(1,(length(commaLocs)+1));
start=1;
for colIdx=1:length(commaLocs)
    fieldNames{colIdx}=tline(start:commaLocs(colIdx)-1);
    start=commaLocs(colIdx)+1;
end
fieldNames{colIdx+1}=tline(start:end);

%Read in the rest of the data (should be numeric)
fieldData=csvread(fileName,1,0);

%Convert the data into a single structure.
[cellW cellH]=size(fieldNames);
numFields=max([cellW cellH]);
%Needs to be a 1xN or Nx1 cell-array.
if ~iscell(fieldNames)
    disp('Needs to be a cell-array');
    return;
elseif cellW ~=numFields && cellH ~=numFields
    disp('Needs to be 1xN or Nx1');
    return;
end

dataSize=size(fieldData);
arrayDim=find(dataSize==numFields);
if isempty(arrayDim)
    disp('Dimensions are wrong');
end

outputStuff=[];
if arrayDim==2
    for idx=1:numFields
        outputStuff.(fieldNames{idx})=fieldData(:,idx);
    end
else
    for idx=1:numFields
        outputStuff.(fieldNames{idx})=fieldData(idx,:);
    end
end
