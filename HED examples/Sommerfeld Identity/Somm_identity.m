% This program computes the Sommerfeld identity
% exp(-jkr) = 4*pi int_{-\infty}^{\infty} J_0(k_p p) k_p 
tic
clear;close all;clf
tol = 1e-15; % tolerance of the routine
num = 200; %Size of the arrays
%% Global Parameters
global i % index number of the distance array
global p % distance
global a % Breakpoint location
global kmax % Switch for TE/TM case (alpha = 0 -> TE, else -> TM)
global maxlev
f = 10e9;
omega = 2*pi*f;
ep1 = 1;
ep2 = 10 - 1i*18;
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;

k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);

a = 2*k1; % Set breakpoint
p = linspace(1e-2,1e2, 200); % Define distance array
% p = linspace(sqrt(2), sqrt(200), 100);
for i = 1 : length(p)
        q = pi/p(i);
        if p(i) < 1
            maxlev = 5;
        else
            maxlev = 15;
        end
        val_1(i) = TanhSinhQuad(0, k1 + .001i, tol); % Integrate upto k through DE
        val_2(i) = TanhSinhQuad(k1 + .001i, a, tol); % Integrate k upto a through DE
%         if p(i) < 1
%             kmax = 10;
%         else
%             kmax = 20;
%         end
        maxlev = 5;
        val_3(i) = PE_Levin(a, tol, q);
        
        val(i) = val_1(i) + val_2(i) + val_3(i);
end

clf
% Individual Contribution
figure (1)

N = 3; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = loglog(p, abs(val_1), 'linewidth',1.3);
hold on
h2 = loglog(p, abs(val_2), 'linewidth',1.3);
h3 = loglog(p, abs(val_3), 'linewidth',1.3);

loglog(p, abs(val_1), 's', 'markersize',4);
loglog(p, abs(val_2), 's', 'markersize',4);
loglog(p, abs(val_3), 's', 'markersize',4);
xlabel('$k_1\rho$','interpreter','latex')
ylabel('$I(z, \rho, \tau)$','interpreter','latex')
legend([h1 h2 h3],{'DE Rule 0 to k', 'DE Rule k to a', 'PE Rule'},...
     'location','northwest');
box on
set(gcf,'color','white');
hold off


% Individual Contribution
figure (2)
f = exp(-1j*k1*p)./(4*pi*p);
N = 2; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = loglog(p, abs(val), 'linewidth',1.3);
hold on
h2 = loglog(p, abs(f), 'linewidth',1.3);
title('Computed Sommerfeld Identity')
loglog(p, abs(val), 's', 'markersize',4);
loglog(p, abs(f), 's', 'markersize',4);
xlabel('$\rho$','interpreter','latex')
ylabel('$\mathcal{S}_0 \left\{  \frac{1}{2jk_{z1}} \right\}$','interpreter','latex')
box on
legend([h1 h2],{'Numerical', 'Analytical'},...
     'location','northeast');
set(gcf,'color','white');
hold off
toc


% Individual Contribution
figure (3)

N = 2; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = loglog(p, abs(val - f), 'linewidth',1.3);
hold on

loglog(p, abs(val - f), 's', 'markersize',4);
xlabel('$\rho$','interpreter','latex')
ylabel('$\vert \mathcal{S}_0 \left\{  \frac{1}{2jk_{z1}} \right\} - \frac{e^{jk_1 \rho}}{4 \pi \rho} \vert$',...
    'interpreter','latex', 'fontsize',11)
title('Difference between analytical and computed values of Sommerfeld Identity')
box on
set(gcf,'color','white');
hold off
toc