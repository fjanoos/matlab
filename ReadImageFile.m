function [I,ijk2ras, file_type] = ReadImageFile( filename )
% dynamically determine the filetype and read it using the appropriate reader
% for .nii files -SPM should be in path !
% for .nrrd file ../nrrd should be in path

    [pathstr, name, ext] = fileparts(filename) ;

    switch ext
        case {'.dcm', '.%03d'}
            [I,ijk2ras]=ReadVolPC16(filename);   
            file_type = 'dicom';
        case '.nii'
            vv = spm_vol(filename);
            ijk2ras = vv.mat; % maybe invert ?
            I = vv.private.dat(:,:,:);
            file_type = 'nifti';
        case '.mhd'
            [I,ijk2ras]=ReadMHDFile(fullfile(pathstr, name));
            file_type = 'metaheader';
        case '.nrrd' % NRRD file     
            vv = nrrdLoadWithMetadata(filename);
            % debug NRRD ijk2ras system
            ijk2ras = [vv.spacedirections,  vv.spaceorigin; 0 0 0 1 ];
            I = vv.data;
            file_type = 'nrrd';
            
            % space field code
            % 1 -> RAS
            % 2 -> LAS
            % 3 -> LPS
            switch vv.space
                case 1
                    display([filename, ' NRRD file in RAS system']);
                case 2
                    display([filename, ' NRRD file in LAS system']);
                case 3
                    display([filename, ' NRRD file in LPS system -changing to RAS']);
                    %ijk2ras = [-1 0 0 0 ; 0 -1 0 0 ; 0 0 1 0; 0 0 0 1]*ijk2ras;
% %                     ijk2ras( 2,3) = - ijk2ras( 2,3);
% %                     ijk2ras( 3,2) = - ijk2ras( 3,2);
            end
    end