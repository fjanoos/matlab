function [n_pts, n_polys] = convertOffToVTK(off_fname, vtk_fname)
% % convert a .off file to a .vtk format for display

    fid = fopen (off_fname);
    
    if fid == -1
        error ('could not open off_fname');
    end
    
    %read in the first line
    %C = textscan( fid, '%s %d %d 0', 1);
    C = textscan( fid, '%s', 1);
    if ~strcmp(C{1} , 'OFF')  
        return
    end
    C = textscan( fid, '%d %d 0', 1);
    n_pts = C{1};
    n_polys = C{2};

    disp(sprintf('off file with %d points and %d polys',n_pts, n_polys));
    
    %read point data
    C = textscan(fid, '%f', n_pts*3);
    pts = []; p = 1;
    for j = 1 :3: 3*n_pts
        pts(p,:) = [C{1}(j),C{1}(j+1),C{1}(j+2)];        
        p = p + 1;
    end

    %read the poly data
    C = textscan(fid, '%d'); %read till end of file
    k = 1; polgon_array = {};
    
    for p = 1 : n_polys
        n_vtx = C{1}(k);
        polygon_array{p}(1) = n_vtx;
        
       for vtx = 1 : n_vtx
           polygon_array{p}(1+vtx) = C{1}(k+vtx);
       end
       k = k+n_vtx+1;
    end

    poly_size = k-1;
    
    fclose( fid);
    
    fid_vtk = fopen (vtk_fname, 'w');
        
    if fid_vtk == -1
        error ('could not open vtk_fname');
    end
    
    %header
    fprintf( fid_vtk, '# vtk DataFile Version 3.0\n');
    fprintf( fid_vtk, 'Generated by le dieu FJ \n');
    fprintf( fid_vtk, 'ASCII\n');
    fprintf( fid_vtk, 'DATASET POLYDATA\n');
    fprintf( fid_vtk, 'POINTS %d float\n', n_pts);
    %write point data
    % [FJ NOTE] - switch the x and y coordinates to match RegisterProstate
    % (buildSurfaceFromLabelMap.m,217)
    for p = 1 : n_pts
        fprintf( fid_vtk, '%f %f %f \n', pts(p,[2 1 3]) );
    end
    fprintf( fid_vtk, 'POLYGONS %d %d\n', n_polys, poly_size);
    %polygon data
    for p = 1 : n_polys
        n_vtx = polygon_array{p}(1);
        fprintf( fid_vtk, '%d ', n_vtx );
       for vtx = 1 : n_vtx
           fprintf( fid_vtk, '%d ', polygon_array{p}(1+vtx) );
       end
       fprintf( fid_vtk, '\n');
    end
    fprintf( fid_vtk, '\nCELL_DATA %d \n', n_polys);
    fprintf( fid_vtk, 'POINT_DATA %d \n', n_pts);
    
    
    fclose( fid);
    
end