% This program computes the Sommerfeld integral for a Vertical Electric
% Dipole
clear all; close all
tic
tol = 1e-6; % tolerance of the routine
num = 60; %Size of the arrays
%% Global Parameters
global i % index number of the distance array
global p % distance
global a % Breakpoint location
global maxlev
f = 10e6;
c = 3e8;
lambda = c/f;
omega = 2*pi*f;

ep1 = 1;
ep2 = 81 - 1i*7190.04;
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;

k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);

a = 2*k1; % Set breakpoint
p = lambda * logspace(0,5, num); % Define distance array
q = pi./p;
% p = 1e3/lambda;
val_1 = zeros(size(p));
val_2 = zeros(size(p));
val_3 = zeros(size(p));
val = zeros(size(p));

for i = 1 : length(p)
        maxlev = 25;
        val_1(i) = TanhSinhQuad(0, k1 + .1i/p(i), tol); % Integrate upto k through DE
        val_2(i) = TanhSinhQuad(k1 + .1i/p(i), a, tol); % Integrate k upto a through DE
        maxlev = 15;
        val_3(i) = PE_Levin(a, tol, q(i)); % Tail through PE Levin with Lucas
    
    val(i) = val_1(i) + val_2(i) + val_3(i);
end
val = -1i*omega/(4*pi)*val; % Ez representation for vertical dipole

% Normalize
val = val/(max(max(abs(val))));
clf
% Individual Contribution
figure (1)

N = 3; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
set(0,'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
h1 = loglog(p/lambda, abs(val_1), 'linewidth',1.3);
hold on
h2 = loglog(p/lambda, abs(val_2), 'linewidth',1.3);
h3 = loglog(p/lambda, abs(val_3), 'linewidth',1.3);

loglog(p/lambda, abs(val_1), 's', 'markersize',4);
loglog(p/lambda, abs(val_2), 's', 'markersize',4);
loglog(p/lambda, abs(val_3), 's', 'markersize',4);

xlabel('$\lambda$','interpreter','latex')
ylabel('Normalized $\vert E_z \vert$','interpreter','latex')
legend([h1 h2 h3],{'DE Rule 0 to k', 'DE Rule k to a', 'PE Rule'},...
     'location','northeast','interpreter','latex');
    title('TM case');
box on
set(gcf,'color','white');
hold off

% Tail Contribution
figure (2)
N = 2; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
set(0,'DefaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
h4 = loglog(p/lambda, abs(val_3), 'linewidth',1.3);
hold on
loglog(p/lambda, abs(val_3), 's', 'markersize',4);
xlabel('$\lambda$','interpreter','latex')
ylabel('Normalized $\vert E_z \vert$','interpreter','latex')

box on
set(gcf,'color','white');
hold off
cleanfigure();

    title('TM case');
    matlab2tikz('filename',sprintf('figures/VED_TM_tail.tex'),'showInfo', false)

clf
% Overall integrals
figure(3)

N = 2; % Number of colors to be used
% Use Brewer-map color scheme 'Set1'
axes('ColorOrder',brewermap(N,'set1'),'NextPlot','replacechildren')
set(0,'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
h5 = loglog(p/lambda, abs(val), 'linewidth',1.3);
hold on
loglog(p/lambda, abs(val), 's', 'markersize',4);
xlabel('$\lambda$','interpreter','latex')
ylabel('Normalized $\vert E_z \vert$','interpreter','latex')
title('TM case');
matlab2tikz('filename',sprintf('figures/VED_TM.tex'),'showInfo', false)
box on
set(gcf,'color','white');
hold off
toc