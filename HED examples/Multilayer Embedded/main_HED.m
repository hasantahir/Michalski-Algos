% This program computes the Sommerfeld integral for a Horizontal Electric
% Dipole
clear all; close all
tic
tol = 1e-9; % tolerance of the routine
num = 200; %Size of the arrays
global maxlev
%% Global Parameters
global i % index number of the distance array
global p % distance
global a % Breakpoint location
global nu % Switch for TE/TM case (alpha = 0 -> TE, else -> TM)
global h
f = .1e12;
omega = 2*pi*f;
lambda = 3e8/f;

% Material Properties
ep1 = 1; % Air
ep2 = 9.7; % GaN layer
ep3 = 11.9; % Silicon base

% EM constants
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;

% Propagations Constants
k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);
k3 = omega*sqrt(mu0*ep0*ep3);


load besselzeros.mat First_15_zeros_J0 First_15_zeros_J1

a = 2*k3; % Set breakpoint
p = linspace(1e-3/k1,1e3/k1, num); % Define distance array
% q = pi./p;
% p = logspace(-2,3, num); % Define distance array
% TE case
nu = 0;
% p = lambda * logspace(-2,3, num); % Define distance array
q = pi./p;

% TM case
% nu = 1;


% Define bessel functions
S_0 = @(kp) besselj(0, kp * p(i));
S_1 = @(kp) besselj(1, kp * p(i));



for i = 1 : length(p)
%     if nu == 0
%         q = First_15_zeros_J0/p(i);
%     else
%         q = First_15_zeros_J1/p(i);
%     end

    % Avoid branch points
    %     if nu == 0 % TE case
    maxlev = 15;
    h = 1;
    val_1(i) = TanhSinhQuad(0, k3 + .0001i, tol); % Integrate upto k through DE
    h = 1;
    val_2(i) = TanhSinhQuad(k3 + .0001i, a, tol); % Integrate k upto a through DE
    h = 1;
    maxlev = 15;
    val_3(i) = PE_Levin(a, tol, q(i)); % Tail through PE Levin with Lucas
    %     else
    %         h = 1;
    %         val_1(i) = TanhSinhQuad(0, k3 + .001i, tol); % Integrate upto k through DE
    %         h = 1;
    %         val_2(i) = TanhSinhQuad(k3 + .001i, a, tol); % Integrate k upto a through DE
    %         h = 1;
    %         val_3(i) = PE_Levin(a, tol, q); % Tail through PE Levin with Lucas
    %     end
    %
    val(i) = val_1(i) + val_2(i) + val_3(i);
end

clf;close all
% Individual Contribution
figure (1)

N = 3; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = loglog(p/lambda, abs(val_1/k1), 'linewidth',1.3);
hold on
h2 = loglog(p/lambda, abs(val_2/k1), 'linewidth',1.3);
h3 = loglog(p/lambda, abs(val_3/k1), 'linewidth',1.3);

loglog(p/lambda, abs(val_1/k1), 's', 'markersize',4);
loglog(p/lambda, abs(val_2/k1), 's', 'markersize',4);
loglog(p/lambda, abs(val_3/k1), 's', 'markersize',4);
xlabel('$k_1\rho$','interpreter','latex')
ylabel('$I(z, \rho, \tau)$','interpreter','latex')
legend([h1 h2 h3],{'DE Rule 0 to k', 'DE Rule k to a', 'PE Rule'},...
    'location','northeast');
if nu == 0
    title('TE case');
else
    title('TM case');
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
cleanfigure();
if nu == 0
    title('TE case');
    matlab2tikz('filename',sprintf('figures/TE_contribution.tex'),'showInfo', false)
else
    title('TM case');
    matlab2tikz('filename',sprintf('figures/TM_contribution.tex'),'showInfo', false)
end



% Tail Contribution
figure (2)
N = 2; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h4 = loglog(p/lambda, abs(val_3/k1), 'linewidth',1.3);
hold on
loglog(p/lambda, abs(val_3/k1), 's', 'markersize',4);
xlabel('$k_1\rho$','interpreter','latex')
ylabel('$I(z, \rho, \tau)$','interpreter','latex')

box on
set(gcf,'color','white');
set(groot,'defaulttextinterpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,...
    'box','on',...
    'FontName','times new roman',...
    'FontSize',12);
hold off
cleanfigure();
if nu == 0
    title('TE case');
    matlab2tikz('filename',sprintf('figures/TE_tail.tex'),'showInfo', false)
else
    title('TM case');
    matlab2tikz('filename',sprintf('figures/TM_tail.tex'),'showInfo', false)
end



% Overall integrals
figure(3)

N = 2; % Number of colors to be used
% Use Brewer-map color scheme 'Set1'
axes('ColorOrder',brewermap(N,'set1'),'NextPlot','replacechildren')

h5 = loglog(p/lambda, abs(val)/k1, 'linewidth',1.3);
hold on
loglog(p/lambda, abs(val)/k1, 's', 'markersize',4);
xlabel('$k_1\rho$','interpreter','latex')
ylabel('$I(z, \rho, \tau)$','interpreter','latex')

if nu == 0
    title('TE case');
else
    title('TM case');
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
toc
cleanfigure();
if nu == 0
    title('TE case');
    matlab2tikz('filename',sprintf('figures/TE.tex'),'showInfo', false)
else
    title('TM case');
    matlab2tikz('filename',sprintf('figures/TM.tex'),'showInfo', false)
end