clear all; close all;clc
tic

f = 1.1e12;
omega = 2*pi*f;
lambda = 3e8/f;

% Material Properties
ep1 = 1; % Air
ep2 = 7.34; % GaN layer
ep3 = 11.9; % Silicon base

% EM constants
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;

% Propagations Constants
k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);
k3 = omega*sqrt(mu0*ep0*ep3);

% Middle Layer thickness
d = lambda;

% Source Location
zp = -lambda/2;

% Layer Heights
z0 = -lambda;
z1 = 0;


% TE/TM switch
nu = 1;

% Wavevector
kz1 = @(kp) sqrt(k1 ^2 - kp .^2);

kz2 = @(kp) sqrt(k2 ^2 - kp .^2);

kz3 = @(kp) sqrt(k3 ^2 - kp .^2);

% Define impedances
if nu == 0
    
    % TE Case
    Z1 = @(kp) (omega*mu0)./kz1(kp);
    Z2 = @(kp) (omega*mu0)./kz2(kp);
    Z3 = @(kp) (omega*mu0)./kz3(kp);
else
    % TM case
    Z1 = @(kp) kz1(kp)./(omega*ep1*ep0);
    Z2 = @(kp) kz2(kp)./(omega*ep2*ep0);
    Z3 = @(kp) kz3(kp)./(omega*ep3*ep0);
end

% Reflection Coefficients
Gamma_left = @(kp)(Z3(kp) - Z2(kp)) ./ (Z3(kp) + Z2(kp)); % Left-looking
Gamma_right = @(kp) (Z1(kp) - Z2(kp)) ./ (Z1(kp) + Z2(kp)); % Right-looking

% Unknown A
A = @(kp) (Gamma_left(kp) .* exp(1i*kz2(kp)*2*z0))./(1 - Gamma_left(kp).*Gamma_right(kp).*exp(-2i * kz2(kp) * d)) ...
    .*( exp(-1i * kz2(kp) * zp) + Gamma_right(kp) .* exp( -1i * kz2(kp) * (2*z1 - zp)));

% Unknown B 
B = @(kp) (Gamma_right(kp) .* exp(-1i*kz2(kp)*2*z1))./(1 - Gamma_left(kp).*Gamma_right(kp).*exp(-2i * kz2(kp) * d)) ...
    .*( exp(+1i * kz2(kp) * zp) + Gamma_left(kp) .* exp( +1i * kz2(kp) * (2*z0 - zp)));

% Denominator
D = @(kp) 1 - Gamma_left(kp).*Gamma_right(kp).*exp(-2i * kz2(kp) * d);
save logical.mat 
%% Pole Finding
% z = AllZeros(D,k2/100,k2*100,1e6); % Pole Locating

% kp Sweep
lxlim = k1/1.5;
uxlim = k3*1.5;
kp = linspace(lxlim,uxlim,1e4);
% 2D kp sweep
[x,y] = meshgrid(linspace(k2/2,k2*2,1e3), linspace(0,1e-2,1e3));
zz = x + 1i*y; % 2d kp
D_lin = D(kp);
D_surf = D(zz);
minn =(min(min(abs(D_lin))));
locate = find( abs(D_lin) == minn);
pole = kp(locate);
bp = k2;
diff = abs( bp - pole);
% clc
% sprintf('Branch Point is at %0.5f',k2)
% sprintf('Pole location is at %0.5f', kp(locate))
% sprintf('Difference of %0.5f', diff)

%% Rational Approximaton of the real and imaginary parts
re_D = real(D_lin);
[rN,rD] = rat(re_D,1e-15);

im_D = imag(D_lin);
[iN,iD] = rat(im_D,1e-15);
% Figures

%% Denominator

figure (1);
N = 4; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% 
h1 = semilogx(kp, real(D_lin), 'linewidth',1.3);
hold on
semilogx(kp, rN./rD, 'linewidth',1.3);

% Imaginary Parts

h1 = semilogx(kp, imag(D_lin), 'linewidth',1.3);
semilogx(kp, iN./iD, 'linewidth',1.3);
% loglog(kp, abs(D), 's', 'markersize',4);
xlabel('$k_{\rho}$','interpreter','latex')
ylabel('$\mathcal{D}$','interpreter','latex')
% legend([h1],{'DE Rule 0 to k'},...
%      'location','northwest');
if nu == 0
    title('Denominator in TE case','interpreter','latex');
else
    title('Denominator in TM case','interpreter','latex');
end
% Decorations
box on
set(gcf,'color','white');
set(groot,'defaulttextinterpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,...
    'box','on',...
    'FontName','times new roman',...
    'FontSize',12);
hold off
xlim([lxlim uxlim])
% matlab2tikz('filename',sprintf('Figures/denominator.tex'),...
%     'showInfo', false,'floatFormat','%.3f');


%% Plot Actual Pole Locations
% figure(2)
% N = 2; % Number of colors to be used
% % Use Brewer-map color scheme
% axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% %
% plot(real(z), (imag(z)), 's', 'markersize',4);
% hold on
% plot(real(k2), imag(k2), 'd', 'markersize',6, 'MarkerFaceColor',[0.5,0.5,0.5])
% %
% xlabel('$\textrm{Real Part of Pole}$','interpreter','latex')
% ylabel('$\textrm{Imaginary Part of Pole}$','interpreter','latex')
% 
% % Decorations
% 
% box on
% set(gcf,'color','white');
% set(groot,'defaulttextinterpreter','latex');
% set(gca,'TickLabelInterpreter', 'latex');
% set(gca,...
%     'box','on',...
%     'FontName','times new roman',...
%     'FontSize',12);
% hold off
% matlab2tikz('filename',sprintf('Figures/poles.tex'),...
%     'showInfo', false,'floatFormat','%.3f');


%% Plot Relative Pole Locations from the Branch Point
% figure(3)
% N = 2; % Number of colors to be used
% % Use Brewer-map color scheme
% axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% %
% plot(real(z) - k2, (imag(z)), 's', 'markersize',4);
% hold on
% plot(real(k2) - k2 , imag(k2), 'd', 'markersize',6,'MarkerFaceColor',[0.5,0.5,0.5]);
% %
% xlabel('$\textrm{Real Relative Distance from Branch Point}$','interpreter','latex')
% ylabel('$\textrm{Imaginary Relative Distance from Branch Point}$','interpreter','latex')
% 
% % Decorations
% 
% box on
% set(gcf,'color','white');
% set(groot,'defaulttextinterpreter','latex');
% set(gca,'TickLabelInterpreter', 'latex');
% set(gca,...
%     'box','on',...
%     'FontName','times new roman',...
%     'FontSize',12);
% hold off
% matlab2tikz('filename',sprintf('Figures/poles_rel_bp.tex'),...
%     'showInfo', false,'floatFormat','%.3f');


%% Surface Plot
% fig = figure(4);
% %
% h = surf(x,y,abs(D_surf),'LineStyle', 'none',...
%     'FaceColor', 'interp');  
% 
% light('Position',[10 0 10],'Style','infinite')
% % Set z-axis to logscale
% set(gca,'Xscale','Log','Yscale','Log','Zscale','Log')
% 
% % Decorations
% 
% % Use Brewermap color schemes
% colormap(brewermap([],'Spectral')) %RdYlGn 'RdYlBu' Spectral
% % Create some lighting for nice figures
% camlight left
% camlight right
% camlight('headlight')
% lighting gouraud
% material metal
% shading interp
% set(gca,...
%     'box','on',...
%     'FontName','times new roman',...
%     'FontSize',12);
% grid on
% set(gcf, 'color','white');
% view([45 30])
% xlim([k2/2 k2*2])
% % Labeling
% xlabel('$\Re {k_{\rho}}$','interpreter','latex')
% ylabel('$\Im {k_{\rho}}$','interpreter','latex')
% % print(fig,'-dpng','figures/surf_denominator','-r600');
%% Denominator Error in rationalization

figure (4);
N = 2; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% 
h1 = loglog(kp, abs(real(D_lin) - rN./rD), 'linewidth',1.3);


% Imaginary Parts
hold on
h1 = loglog(kp, abs(imag(D_lin) - iN./iD), 'linewidth',1.3);

% loglog(kp, abs(D), 's', 'markersize',4);
xlabel('$k_{\rho}$','interpreter','latex')
ylabel('$\mathcal{D}$','interpreter','latex')
% legend([h1],{'DE Rule 0 to k'},...
%      'location','northwest');
if nu == 0
    title('Denominator Error in TE case','interpreter','latex');
else
    title('Denominator Error in TM case','interpreter','latex');
end
% Decorations
box on
set(gcf,'color','white');
set(groot,'defaulttextinterpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,...
    'box','on',...
    'FontName','times new roman',...
    'FontSize',12);
hold off
xlim([lxlim uxlim])
toc


figure (5)

N = 4; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% 
semilogx(kp, real(Gamma_right(kp)), 'linewidth',1.3);
hold on
semilogx(kp, real(Gamma_left(kp)), 'linewidth',1.3);
semilogx(kp, imag(Gamma_left(kp)), 'linewidth',1.3)
semilogx(kp, imag(Gamma_right(kp)), 'linewidth',1.3)

xlabel('$k_{\rho}$','interpreter','latex')
ylabel('$\mathcal{D}$','interpreter','latex')

if nu == 0
    title('$\Gamma$ functions in TE case','interpreter','latex');
else
    title('$\Gamma$ functions in TM case','interpreter','latex');
end
% Decorations
box on
set(gcf,'color','white');
set(groot,'defaulttextinterpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,...
    'box','on',...
    'FontName','times new roman',...
    'FontSize',12);
hold off
xlim([lxlim uxlim])