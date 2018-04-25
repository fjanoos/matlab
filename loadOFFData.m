function [pts, polygon_array] = loadOFFData(off_fname)
% read the polygon mesh from the off file

fid = fopen (off_fname);
    
%read in the first line
C = textscan( fid, '%s %d %d 0', 1);

if strcmp(C{1}{1} , 'OFF')  
    n_pts = C{2}(1);
    n_polys = C{3}(1);

    disp(sprintf('off file with %d points and %d polys',n_pts, n_polys));
else
    return
end

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

end