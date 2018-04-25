function WriteMHDFile(filebase, I, transform_matrix, elem_size, offset)
% writes an metaheader type file
% Parameters:
%   I - 3D Image 
%   transform_matrix - either  a 4x4 matrix
%                      if it is a 3x3 transformation matrix, then need to
%                      specify the next two parameters
%   elem_size        - a 3x1 scaling vector
%   offset           - a 3x1 offset vector



disp(sprintf('Writing %s',filebase));

if length(size(I)) ~= 3
    error 'support only for 3d images as of now';
end

hdr_fname = [filebase, '.mhd'];
vol_fname = [filebase, '.raw'];
hdr_fid = fopen( hdr_fname, 'w');
vol_fid = fopen( vol_fname, 'w');

[pth,fbase,ext] = fileparts(filebase);


if ~exist(hdr_fname, 'file') || ~exist(vol_fname, 'file')
    error( 'could not create files ! ');
    return;
end


if size(transform_matrix) == [4,4]
    %error('4x4 feature not supported');
    transform_matrix = transform_matrix/transform_matrix(4,4); %homogenize
    offset= transform_matrix(1:3,4);
    if any( transform_matrix(4,1:3) ~= 0 )
        error ('error in tranform matrix' );
    end
    transform_matrix(:,4) = [];
    transform_matrix(4,:) = [];    
    elem_size = sqrt(diag( transform_matrix'*transform_matrix ));
    transform_matrix = transform_matrix/diag(elem_size);
end
   
fprintf(hdr_fid, ['ObjectType = Image\n',...
                  'NDims = 3\n',...
                  'BinaryData = True\n',...
                    'BinaryDataByteOrderMSB = False\n',...
                'CompressedData = False\n']);
fprintf(hdr_fid, 'TransformMatrix = %g %g %g %g %g %g %g %g %g\n', ...
                    [transform_matrix(:,1),transform_matrix(:,2),transform_matrix(:,3)]);
fprintf(hdr_fid, 'Offset = %g %g %g\n', offset);
fprintf(hdr_fid, 'CenterOfRotation = 0 0 0\n');
fprintf(hdr_fid, 'AnatomicalOrientation = RAS\n');
fprintf(hdr_fid, 'ElementSpacing = %g %d %d\n', elem_size);
fprintf(hdr_fid, 'DimSize = %g %g %g\n', size(I));
fprintf(hdr_fid, 'ElementSize = %g %g %g\n', elem_size);

switch( class(I) )
    case {'int8', 'int16' , 'int32', 'uint8' , 'uint16','uint32'}
        fprintf(hdr_fid, 'ElementType = MET_SHORT\n');
        fwrite(vol_fid, uint16(I(:)), 'uint16');
        
    case {'single', 'float' , 'float32', 'real*4', 'double'}
        fprintf(hdr_fid, 'ElementType = MET_FLOAT\n');
        fwrite(vol_fid, single(I(:)), 'float32');
        
    otherwise
        error('unsupported datatype');
end
        

fprintf(hdr_fid, 'ElementDataFile = %s.raw', fbase);

fclose(hdr_fid);
fclose(vol_fid);

