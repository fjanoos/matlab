function h = fj_subplot( m,n, idx )
% a tighter plot than matlabs
scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/4])

h_spacing = 0.1/n;
v_spacing = 0.1/m;
height = 1/m ;
width = 1/n;

r = floor((idx-1) / n)% + 1
c = mod( idx-1 , n )% 

h = axes( 'Position', [c*width+h_spacing r*height+v_spacing...
                            width-h_spacing height-2*v_spacing ] );
                       
drawnow;