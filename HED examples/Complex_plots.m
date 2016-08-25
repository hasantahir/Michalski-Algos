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

S_0 = @(kp) G_1(kp) .* besselj(0, kp * 1e1) .* kp;
S_1 = @(kp) G_2(kp) .* besselj(1, kp * 1e1) .* kp;

[X, Y] = meshgrid(0:0.008:2,   0:0.005:.1);
z = X + 1i*Y;

% Sommerfeld Integral kernels
f_0 = S_0(z);
f_1 = S_1(z);

% Surface plot
surf(X,Y,abs(f_1),'LineStyle', 'none', 'FaceColor', 'interp');   

% Set z-axis to logscale
set(gca,'Xscale','Lin','Yscale','Lin','Zscale','Log')

% Use Brewermap color schemes
colormap(brewermap([],'RdYlBu'))

% Create some lighting for nice figures
camlight left
camlight right
camlight('headlight')

% Labeling
xlabel('\Re \rho')
ylabel('\Im \rho')
zlabel('$\mathcal{S}_0$', 'interpreter', 'latex')
shading interp
