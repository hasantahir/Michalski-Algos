close all
figure(1)
N = 4; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren');
Colord = get(gca, 'ColorOrder');
%
plot(real(root)/k1 , (imag(root)/k1), 's', 'markersize',3.5,...
    'MarkerFaceColor',Colord(1,:));
hold on
plot(real(lossless)/k1 , (imag(lossless)/k1), 'o', 'markersize',3.5,...
    'MarkerFaceColor',Colord(2,:));
plot(real(k1)/k1  , imag(k1)/k1 , '<', 'markersize',7,...
    'MarkerFaceColor',Colord(3,:));

plot(real(k2)/k1  , imag(k2)/k1 , 'p', 'markersize',7,...
    'MarkerFaceColor',Colord(4,:));
%
xlabel('$\Re(\textrm{k}_{\rho})$','interpreter','latex')
ylabel('$\Im(\textrm{k}_{\rho})$','interpreter','latex')
legend('With PEC base','Without PEC base','k_{air}','k_{GaN}',...
    'Location','northeast','Orientation','horizontal');
if nu == 0
    title('Zeros of denominator $\mathcal D$');
else
    title('Zeros of denominator $\mathcal D$');
end
box on
set(gcf,'color','white');
set(groot,'defaulttextinterpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,...
    'box','on',...
    'FontName','times new roman',...
    'FontSize',12);
hold off
ylim([-10 10])
% Save tikz figure
% cleanfigure();
if nu == 0
    matlab2tikz('filename',sprintf('figures/TE_pole_comparison_d_%d.tex',floor(d*1e7)),'showInfo', false);
else
    matlab2tikz('filename',sprintf('figures/TM_pole_comparison_d_%d.tex',floor(d*1e7)),'showInfo', false);
end
