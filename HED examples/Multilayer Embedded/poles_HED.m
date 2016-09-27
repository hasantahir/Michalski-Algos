% This program computes the poles for a Horizontal Electric
% Dipole embedded in multilayer environment
clear; close all
tic



%% Global Parameters
num = 1e3; %Size of the arrays
global nu % Switch for TE/TM case (alpha = 0 -> TE, else -> TM)
global f
global ep1 ep2 ep3
global k1 k2 k3
global mu0 ep0 


f = 1.1e12;
omega = 2*pi*f;
lambda = 3e8/f;

% Material Properties
ep1 = 1; % Air
ep2 = 9.7 + .001i;%- .01i; % GaN layer
ep3 = 11.9; % Silicon base

% EM constants
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;

% Propagations Constants
k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);
k3 = omega*sqrt(mu0*ep0*ep3);

lxlim = k1/1.50;
uxlim = k3*1.50;
p = linspace(lxlim,uxlim,num);
kp = meshgrid(p,p);

% TE case
% nu = 0;

% TM case
nu = 1;

for i = 1 : length(p)

D(i) = Denom(p(i));

end
% for i = 1 : length(p)
%     for j = 1 :length(p)
%         D(i,j) = Denom(kp(i,j));
%     end
% end
D = D';

% Figures

%% Denominator

figure (1);
N = 2; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% 
h1 = semilogx(p, real(D), 'linewidth',1.3);
hold on


% Imaginary Parts

semilogx(p, imag(D), 'linewidth',1.3);

% loglog(kp, abs(D), 's', 'markersize',4);
xlabel('$k_{\rho}$','interpreter','latex')
ylabel('$\mathcal{D}$','interpreter','latex')
legend('Real Part','Imaginary Part','location','northeast')
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


toc