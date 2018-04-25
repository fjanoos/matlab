function [n_pts, n_polys] = convertVTKToOff(vtk_fname, off_fname)
%% convert a vtk ascii file to .off format
display( ['Running ', mfilename ]);

fid = fopen (vtk_fname);

% read in the file version id # vtk DataFile Version 3.0

C = textscan( fid, '%s', 5);

if strcmp(C{1}{1}, '#') & strcmp(C{1}{2}, 'vtk') & ...
        strcmp(C{1}{3} ,'DataFile') &strcmp(C{1}{4}, 'Version')
    display( ['vtk data file version ' C{1}{5}] );
end

%read the next line
C = textscan( fid, '%s', 2,  'delimiter', '\n');

str = '';
for k = 1:length(C)
   str = C{1}{k};
   disp(str);
end

%read the next 2 line ASCII && 'DATASET POLYDATA'
C = textscan( fid, '%s', 2,  'delimiter', '\n');

if strcmp (C{1}{1} , 'ASCII' ) & strcmp (C{1}{2} , 'DATASET POLYDATA' )
    disp('ASCII polydata file')
else
    return
end


%read the next line POINTS n_pts datatype
C = textscan( fid, '%s %d %s', 1);

if strcmp(C{1}{1} , 'POINTS')  
    n_pts = C{2};
    
    disp(sprintf('number points %d of %s type',n_pts, C{3}{1}));
else
    return
end

% read in the points
C = textscan(fid, '%f', n_pts*3);

pts = []; j = 1;
for k = 1:3:length(C{1})
   pts(j,:) = [C{1}(k),C{1}(k+1),C{1}(k+2)];
   j = j+1;
end
%disp(pts);


%read the next line POLYGONS ncells size
C = textscan( fid, '%s %d %d', 1);

if strcmp(C{1}{1} , 'POLYGONS')  
    n_polys = C{2}(1);
    poly_list_size = C{3}(1);
    %disp(sprintf('number polys %d, size %d',n_polys, poly_list_size));
else
    return
end


% read in the polygons
C = textscan(fid, '%d', poly_list_size);

polgon_array = {}; j = 1; k = 1;
while k < length(C{1})
   n_vtx = C{1}(k);
   polygon_array{j}(1) = n_vtx;
   for vtx = 1 : n_vtx
       polygon_array{j}(1+vtx) = C{1}(k+vtx);
   end
 %  disp(polygon_array{j});
   j = j+1;
   k = k+n_vtx+1;
end

%ignore everything else
% % %read in the next set of strings until i hit normal
% % while(1)
% %     C = textscan( fid, '%s', 1);
% %     if ( strcmp(C{1}{1}, 'NORMALS' ) )
% %         %read till end of line
% %         C = textscan( fid, '%s', 1,  'delimiter', '\n');
% %         %disp(C{1});
% %         % read in the normals (for each point)
% %         C = textscan(fid, '%f', n_pts*3);
% %         normals = []; j = 1;
% %         for k = 1:3:length(C{1})
% %            normals(j,:) = [C{1}(k),C{1}(k+1),C{1}(k+2)];
% %            j = j+1;
% %         end
% %        % disp(normals);
% %         break;
% %     end
% % end

fclose(fid);

fid_off = fopen (off_fname, 'w');

%header
fprintf( fid_off, 'OFF %d %d 0\n', n_pts, n_polys);
%point data
for p = 1 : n_pts
    fprintf( fid_off, '%f %f %f \n', pts(p,:) );
end
%polygon data
for p = 1 : n_polys
    n_vtx = polygon_array{p}(1);
    fprintf( fid_off, '%d ', n_vtx );
   for vtx = 1 : n_vtx
       fprintf( fid_off, '%d ', polygon_array{p}(1+vtx) );
   end
   fprintf( fid_off, '\n');
end

fclose(fid_off);


    
end