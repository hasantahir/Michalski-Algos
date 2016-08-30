% This program computes the Sommerfeld identity
% exp(-jkr) = 4*pi int_{-\infty}^{\infty} J_0(k_p p) k_p 
tic
clear
tol = 1e-15; % tolerance of the routine
num = 200; %Size of the arrays
%% Global Parameters
global i % index number of the distance array
global p % distance
global a % Breakpoint location
global nu % Switch for TE/TM case (alpha = 0 -> TE, else -> TM)
f = 10e9;
omega = 2*pi*f;
ep1 = 1;
ep2 = 10 - 1i*18;
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;

k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);

load besselzeros.mat First_15_zeros_J0 First_15_zeros_J1

nu = 0;
a = 2*k1; % Set breakpoint
p = linspace(1e-2,1e0, 100); % Define distance array
% p = linspace(sqrt(2), sqrt(200), 100);
for i = 1 : length(p)
%         q = First_15_zeros_J0/p(100);
        q = pi/p(i);
        val_1(i) = TanhSinhQuad(0, k1 + .01i, tol); % Integrate upto k through DE
        val_2(i) = TanhSinhQuad(k1 + .01i, a, tol); % Integrate k upto a through DE
        val_3(i) = 0;%TanhSinhQuad(a, q(1), tol); % Integrate k upto a through DE
        val_4(i) = PE_Levin(a, tol, q);
        
        val(i) = val_1(i) + val_2(i) + val_3(i) + val_4(i);
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
% legend([h1 h2 h3],{'DE Rule 0 to k', 'DE Rule k to a', 'PE Rule'},...
%      'location','northwest');
if nu == 0
    title('TE case');
else
    title('TM case');
end
box on
set(gcf,'color','white');
hold off


clf
% Individual Contribution
figure (2)

N = 2; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = loglog(p, abs(val), 'linewidth',1.3);
hold on

loglog(p, abs(val), 's', 'markersize',4);
xlabel('$k_1\rho$','interpreter','latex')
ylabel('$I(z, \rho, \tau)$','interpreter','latex')
% legend([h1 h2 h3],{'DE Rule 0 to k', 'DE Rule k to a', 'PE Rule'},...
%      'location','northwest');
if nu == 0
    title('TE case');
else
    title('TM case');
end
box on
set(gcf,'color','white');
hold off
toc