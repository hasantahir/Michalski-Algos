x = linspace(.95*k_air,1.15*k_air,300);
y = linspace(-.0005*k_air,.0005*k_air,300);
[x,y] = meshgrid(x,y);
kp = x + 1i*y;
% D = @(kp) 54.0 + kp.*(44.0 + kp.*(20.0 - kp.*(3.0-kp)));
% D = @(kp)kp.^50 + kp.^12 - 5*sin(20*kp).*cos(12*kp) - 1;
% DD = D(kp);
[D,dD] = FZ(kp);
h(1) = surf(x,y,abs(D),'FaceColor','interp','EdgeColor','interp');
colormap(brewermap([],'RdYlBu')) %RdYlGn 'RdYlBu' Spectral
freezeColors %freeze this plot's colormap
% hold on
% h(2) = surf(x,y,imag(D),'FaceColor','interp','EdgeColor','interp');
% colormap(brewermap([],'Spectral')) %RdYlGn 'RdYlBu' Spectral
% freezeColors %freeze this plot's colormap
% % h(3) = surf(x,y,zeros(size(kp)),'FaceColor','interp','EdgeColor','interp');
% colormap gray
% hold off
% view([0 90])
% % xlim([-5 5])
% % ylim([-5 5])