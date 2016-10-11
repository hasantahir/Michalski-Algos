close all;clf
clear;
f = 10e9;
omega = 2*pi*f;
ep1 = 1;
ep2 = 10 - 1i*18;
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;
k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);

kz1 = @(kp) sqrt(k1 ^2 - kp .^2);
kz2 = @(kp) sqrt(k2 ^2 - kp .^2);

gamma_1h = @(kp) -(kz2(kp) - kz1(kp)) ./ (kz2(kp) + kz1(kp));
gamma_1e = @(kp) (kz2(kp) / ep2 - kz1(kp)) ./ (kz2(kp) / ep2 + kz1(kp));

G_1 = @(kp) k1 * (gamma_1h(kp)) ./ (1i * kz1(kp));
G_2 = @(kp) k1 ./ (kp .* (gamma_1e(kp) - gamma_1h(kp)));

S_1 = @(kp) G_1(kp) .* besselj(0, kp * 1e1) .* kp;
S_2 = @(kp) G_2(kp) .* besselj(1, kp * 1e1) .* kp;

% Bessel Function
J_0 = @(kp) besselj(0, kp * 1e1) .* kp;
J_1 = @(kp) besselj(1, kp * 1e1) .* kp;

% Bessel Function
[X, Y] = meshgrid(linspace(0,2,500),   linspace(0,1,500));

% [X, Y] = meshgrid(0:0.008:2,   -1:0.005:1);

% For Spectral Greens function
% [X, Y] = meshgrid(linspace(208,212,500), linspace(0,.5,500));
% X = single(X);
% Y = single(Y);
z = X + 1i*Y;

% Sommerfeld Integral kernels
f_0 = kz1(z);
f_1 = kz2(z);

f_2 = gamma_1h(z);
f_3 = gamma_1e(z);

f_4 = G_1(z);
f_5 = G_2(z);

f_6 = S_1(z);
f_7 = S_2(z);

f_8 = J_0(z);
f_9 = J_1(z);

plt = abs(f_9);
fig = figure;
%% Surface plot
h = surf(X,Y,plt,'LineStyle', 'none',...
    'FaceColor', 'interp');  
%%

plt_min = min(min(abs(plt)));
plt_max = max(max(abs(plt)));
light('Position',[10 0 10],'Style','infinite')
% Set z-axis to logscale
set(gca,'Xscale','Lin','Yscale','Log','Zscale','Log')

% Use Brewermap color schemes
colormap(brewermap([],'Spectral')) %RdYlGn 'RdYlBu' Spectral
% Create some lighting for nice figures
camlight left
camlight right
camlight('headlight')
lighting gouraud
material metal
shading interp
set(gca,...
    'box','on',...
    'FontName','times new roman',...
    'FontSize',12);
grid on
set(gcf, 'color','white');
view([45 30])
ylim([1e-2 1])
% Labeling
xlabel('$\Re {k_{\rho}}$','interpreter','latex')
ylabel('$\Im {k_{\rho}}$','interpreter','latex')
% caxis([abs((plt_min)) .04*abs(log(plt_max))])


if plt == abs(f_8)
%     title('$J_0(z)$', 'interpreter', 'latex');
%     matlab2tikz('filename',sprintf('figures/J_0_comp.tex'),...
%         'showInfo', false,'floatFormat','%.3f');
    print(fig,'-dpng','figures/J_0_comp','-r600')
elseif plt == abs(f_9)
%     title('$J_1(z)$', 'interpreter', 'latex');
%     matlab2tikz('filename',sprintf('figures/J_1_comp.tex'),...
%         'showInfo', false,'floatFormat','%.3f');
    print(fig,'-dpng','figures/J_1_comp','-r600')
elseif plt == abs(f_6)
    print(fig,'-dpng','figures/S_1_comp','-r600')
end