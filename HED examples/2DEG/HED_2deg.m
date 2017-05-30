% This program computes the Sommerfeld integral for a Horizontal Electric
% Dipole
clear all; close all
tic
tol = 1e-15; % tolerance of the routine
num = 60; %Size of the arrays
%% Global Parameters
global i % index number of the distance array
global p % distance
global a % Breakpoint location
global nu % Switch for TE/TM case (alpha = 0 -> TE, else -> TM)
global h
% f = 10e9;
% omega = 2*pi*f;
% ep1 = 1;
% ep2 = 10 - 1i*18;
% mu0 = 4*pi*1e-7;
% ep0 = 8.854e-12;
% 
% k1 = omega*sqrt(mu0*ep0*ep1);
% k2 = omega*sqrt(mu0*ep0*ep2);

% STO 
% f = 2.8e12;
% omega = 2*pi*f;
% ep1 = 1;
% ep2 = -90 + 1i*577;
% mu0 = 4*pi*1e-7;
% ep0 = 8.854e-12;
% 
% k1 = omega*sqrt(mu0*ep0*ep1);
% k2 = omega*sqrt(mu0*ep0*ep2);

% GaAs
% f = 8.4e12;
% omega = 2*pi*f;
% ep1 = 1;
% ep2 = -12.5 + 1i*3.2;
% mu0 = 4*pi*1e-7;
% ep0 = 8.854e-12;

% Sea water
f = 10e6;
c = 3e8;
omega = 2*pi*f;
lambda = c/f;
ep1 = 1;
ep2 = 81 -1i*8192;
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;

k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);

load bess_zeros.mat j0 j1

a = 2*k1; % Set breakpoint
p = linspace(1e-0/k1,1e4/k1, num); % Define distance array
% p = lambda * logspace(0,4, num); % Define distance array
q = pi./p;
% TE case
nu = 0;


% TM case
% nu = 1;


% Define bessel functions
S_0 = @(kp) besselj(0, kp * p(i));
S_1 = @(kp) besselj(1, kp * p(i));



for i = 1 : length(p)
    if nu == 0
        q = j0/p(i);
    else
        q = j1/p(i);
    end
    % Avoid branch points
    if nu == 0 % TE case
        
        h = 1;
        val_1(i) = TanhSinhQuad(0, a, tol); % Integrate upto k through DE
%         h = 1;
%         val_2(i) = TanhSinhQuad(k1, a, tol); % Integrate k upto a through DE
        h = 1;
        val_3(i) = PE_Levin(a, tol, q); % Tail through PE Levin with Lucas
    else
        h = 1;
        val_1(i) = TanhSinhQuad(0, a, tol); % Integrate upto k through DE
%         h = 1;
%         val_2(i) = TanhSinhQuad(k1, a, tol); % Integrate k upto a through DE
        h = 1;
        val_3(i) = PE_Levin(a, tol, q); % Tail through PE Levin with Lucas
    end
    
    val(i) = val_1(i) + val_3(i);
end

clf
% Individual Contribution
% figure (1)
% 
% N = 3; % Number of colors to be used
% % Use Brewer-map color scheme
% axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% % h1 = loglog(p*k1, abs(val_1/k1), 'linewidth',1.3);
% hold on
% % h2 = loglog(p*k1, abs(val_2/k1), 'linewidth',1.3);
% h3 = loglog(p*k1, abs(val_3/k1), 'linewidth',1.3);
% 
% % loglog(p*k1, abs(val_1/k1), 's', 'markersize',4);
% % loglog(p*k1, abs(val_2/k1), 's', 'markersize',4);
% loglog(p*k1, abs(val_3/k1), 's', 'markersize',4);
% xlabel('$\rho$','interpreter','latex')
% ylabel('$I(z, \rho, \tau)$','interpreter','latex')
% % legend([h1 h2 h3],{'DE Rule 0 to k', 'DE Rule k to a', 'PE Rule'},...
% %      'location','northwest');
% if nu == 0
%     title('TE case');
% else
%     title('TM case');
% end
% box on
% set(gcf,'color','white');
% hold off

% Tail Contribution
figure (2)
N = 2; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h4 = loglog(p*k1, abs(val_3/k1), 'linewidth',1.3);
hold on
loglog(p*k1, abs(val_3/k1), 's', 'markersize',4);
xlabel('$k_1\rho$','interpreter','latex')
ylabel('$I(z, \rho, \tau)$','interpreter','latex')

box on
set(gcf,'color','white');
hold off
% cleanfigure();
if nu == 0
    ylabel('$E_{\pho}$','interpreter','latex');
    matlab2tikz('filename',sprintf('figures/E_p.tex'),'showInfo', false)
else
    ylabel('$E_{z}$','interpreter','latex');
    matlab2tikz('filename',sprintf('figures/E_z.tex'),'showInfo', false)
end



% % Overall integrals
figure(3)

N = 2; % Number of colors to be used
% Use Brewer-map color scheme 'Set1'
axes('ColorOrder',brewermap(N,'set1'),'NextPlot','replacechildren')

h5 = loglog(p*k1, abs(val)/k1, 'linewidth',1.3);
hold on
loglog(p*k1, abs(val)/k1, 's', 'markersize',4);
xlabel('$k_1\rho$','interpreter','latex')
ylabel('$I(z, \rho, \tau)$','interpreter','latex')

if nu == 0
    title('TE case');
else
    title('TM case');
end
box on
set(gcf,'color','white');
% hold off
toc