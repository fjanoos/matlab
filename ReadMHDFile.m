function [I, ijk2ras]=ReadMHDFile(filespec)
% reads an metaheader type file
% returns an image volume in ijk coords and coord transform matrix
% returns int16 volume
% if [M,N,P]=size(image), then
% i is columns of image 1...N
% j is rows of image 1...M
% k is plane of image 1..P
display( ['Running ', mfilename ]);
disp(sprintf('Reading %s',filespec));

hdr_fname = [filespec, '.mhd'];
vol_fname = [filespec, '.raw'];
if ~exist(hdr_fname, 'file') || ~exist(vol_fname, 'file')
    error( 'could not locate file ! ');
    return;
end

hdr_fid = fopen( hdr_fname, 'r');

ijk2ras = eye(4);
offset = zeros(4,4);
scale = eye(4);
dims = [0,0,0];
dtype = 'MET_FLOAT';
I = [];
rai2ras = eye(4);
% % rai2ras(1,1) = -1; %convert into lps
% % rai2ras(2,2) = -1;
rai2ras(3,3) = 1;
rai2ras_flag = 0;

tline = fgetl(hdr_fid);
while ischar(tline)
    %disp(tline);
    [word, r] = strtok(tline);
    switch word
        case 'NDims'
            [n, r] = strtok(r);
            if strtok(r) ~= '3'
                error( 'incorrect number of dimesions');
                return;
            end 
        case 'TransformMatrix'
            [n, r] = strtok(r);
            for j = 1:3
               for i = 1:3
                [n, r] = strtok(r);
                ijk2ras(i,j) = str2num(n);
               end
            end
        case 'Offset'
            [n, r] = strtok(r);
            for i = 1:3               
                [n, r] = strtok(r);
                offset(i,4) = str2num(n);
            end
        case 'AnatomicalOrientation'
            [n, r] = strtok(r);
            if strtok(r) == 'RAI'
                %error( 'not RAI system - consider adding code for flipping');
                rai2ras_flag = 1;
            end 
         case 'ElementSpacing'
            [n, r] = strtok(r);
            for i = 1:3               
                [n, r] = strtok(r);
                scale(i,i) = str2num(n);
            end
         case 'DimSize'
            [n, r] = strtok(r);
            for i = 1:3               
                [n, r] = strtok(r);
                dims(i) = str2num(n);
            end   
         case 'ElementType'
            [n, r] = strtok(r);
            dtype = strtok(r)
            if strcmpi(strtok(r), 'MET_FLOAT')
               disp('floating datatype');
            end   
           
        otherwise
            %disp(tline);            
             
    end
    
    tline = fgetl(hdr_fid);
end

ijk2ras = rai2ras*(ijk2ras*scale+ offset);

vol_fid = fopen( vol_fname, 'r');
switch (dtype)
    case 'MET_FLOAT'
        vol_data = fread(vol_fid, dims(1)*dims(2)*dims(3), 'float32');
    case 'MET_SHORT'
        vol_data = fread(vol_fid, dims(1)*dims(2)*dims(3), 'short');
end

I = int16(reshape(vol_data, dims));


