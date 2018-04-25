function h1 = fj_plot( wd, ht)
% a tighter plot than matlabs
%unit in points
h=figure('Units','points',...
            'Position',[100, 100, wd ht])
set( h, 'PaperUnits', 'points', 'PaperSize', [wd,ht]);        

% h_spacing = 0.1/n;
% v_spacing = 0.1/m;
% height = 1/m ;
% width = 1/n;
% 
% r = floor((idx-1) / n)% + 1
% c = mod( idx-1 , n )% 
% 
h1 = axes('Units','points', 'Position', [18, 38, wd-24, ht-42] );
%h1 = axes('Units','points', 'Position', [38, 28, wd-74, ht-62] );
set( h1, 'FontSize', 10, 'FontName', 'Times New Roman' ),

drawnow;
