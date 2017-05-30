clear;close all;tic
% f = 1e12;
f = 5e9;
omega = 2*pi*f;
lambda = 3e8/f;
num = 3e0; %Size of the arrays
tol = 1e-12;
% Example Validations

% Material Properties
ep1 = 12; 
ep2 = 4 ; 
ep3 = 2.1;
ep4 = 1;


% EM constants
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;

% Propagations Constants
k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);
k3 = omega*sqrt(mu0*ep0*ep3);
k4 = omega*sqrt(mu0*ep0*ep4);




% Source Location
zp = 5e-3;

% Layer Heights
z1 = 0;
z2 = 4e-3;
z3 = 6e-3;
z4 = 10e-3;
dp = zp - z2;
dpp = z3 - zp;

% Layer thickness
d1 = z2 - z1;
d2 = z3 - z2;
d3 = z4 - z3;

% TE/TM switch
nu = 1;

kz1 = @(kp) sqrt(k1 ^2  - kp .^2);

kz2 = @(kp) sqrt(k2 ^2  - kp .^2);

kz3 = @(kp) sqrt(k3 ^2- kp .^2);

kz4 = @(kp) sqrtbr(k4 ^2- kp .^2);

% Define impedances
if nu == 0
    
    % TE Case
    Z1 = @(kp) omega./kz1(kp);
    Z2 = @(kp) omega./kz2(kp);
    Z3 = @(kp) omega./kz3(kp);
    Z4 = @(kp) omega./kz4(kp);
else
    % TM case
    Z1 = @(kp) kz1(kp)./(omega*ep1*ep0);
    Z2 = @(kp) kz2(kp)./(omega*ep2*ep0);
    Z3 = @(kp) kz3(kp)./(omega*ep3*ep0);
    Z4 = @(kp) kz4(kp)./(omega*ep4*ep0);
end


Gamma_12 = @(kp) (Z1(kp) - Z2(kp))./(Z1(kp) + Z2(kp));
Gamma_32 = @(kp) (Z3(kp) - Z2(kp))./(Z3(kp) + Z2(kp));
Gamma_43 = @(kp) (Z4(kp) - Z3(kp))./(Z4(kp) + Z3(kp));


% % Reflection Coefficients
Gamma_left = @(kp) (Gamma_12(kp) + (-1)*exp(-2i*kz1(kp)*d1))...
    ./(1 + Gamma_12(kp)*(-1).*exp(-2i*kz1(kp)*d1)); % Left-looking
Gamma_right = @(kp) (Gamma_32(kp) + Gamma_43(kp).*exp(-2i*kz3(kp)*d3))...
    ./(1 + Gamma_32(kp).*Gamma_43(kp).*exp(-2i*kz3(kp)*d3)); % Right-looking

% Unknown A
A = @(kp) (Gamma_left(kp) .* exp(1i*kz2(kp)*2*z2))./(1 - Gamma_left(kp).*Gamma_right(kp).*exp(-2i * kz2(kp) * d2)) ...
    .*( exp(-1i * kz2(kp) * zp) + Gamma_right(kp) .* exp( -1i * kz2(kp) * (2*z3 - zp)));

% Unknown B
B = @(kp) (Gamma_right(kp) .* exp(-1i*kz2(kp)*2*z3))./(1 - Gamma_left(kp).*Gamma_right(kp).*exp(-2i * kz2(kp) * d2)) ...
    .*( exp(+1i * kz2(kp) * zp) + Gamma_left(kp) .* exp( +1i * kz2(kp) * (2*z2 - zp)));

% Denominator
D = @(kp) 1 - Gamma_left(kp).*Gamma_right(kp).*exp(-2i * kz2(kp) * d2);

lxlim = 1*k4;
uxlim = 5.5*k4;
p = linspace(lxlim,uxlim,num);
root = [];
for i = 1 : length(p)
    r = newtzero(D,p(i));
    root = vertcat(root,r);
end
% Sort the array
root = sort(root);
% Clean up roots by weeding out too close values
% if ~isempty(root)
%     cnt = 1;  % Counter for while loop.
%     
%     while ~isempty(root)
%         vct = abs(root - root(1)) < 1e1; % Minimum spacing between roots.
%         C = root(vct);  % C has roots grouped close together.
%         [idx,idx] = min(abs(D(C)));  % Pick the best root per group.
%         rt(cnt) = C(idx); %  Most root vectors are small.
%         root(vct) = []; % Deplete the pool of roots.
%         cnt = cnt + 1;  % Increment the counter.
%     end
%     root = sort(rt).';  % return a nice, sorted column vector
% end
roots = root((real(root))>k4);
% roots = roots;
if ~isempty(roots)
    cnt = 1;  % Counter for while loop.
    
    while ~isempty(roots)
        vct = abs(roots - roots(1)) < 1e-4; % Minimum spacing between roots.
        C = roots(vct);  % C has roots grouped close together.
        [idx,idx] = min(abs(D(C)));  % Pick the best root per group.
        rt(cnt) = C(idx); %  Most root vectors are small.
        roots(vct) = []; % Deplete the pool of roots.
        cnt = cnt + 1;  % Increment the counter.
    end
    roots = sort(rt).';  % return a nice, sorted column vector
end
%% Physical Roots
% Real_roots = root(real(root)>k1);
% % Plot
figure(1)
N = 5; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
Colord = get(gca, 'ColorOrder');
plot(real(root)/k4, imag(root)/k4, 's', 'markersize',4,...
    'MarkerFaceColor',Colord(1,:));

hold on
plot(real(k4)/k4, imag(k4)/k4,'d', 'markersize',4,...
    'MarkerFaceColor',Colord(2,:));

% plot(real(k1)/k1 , imag(k1)/k1, 'd', 'markersize',4,...
%     'MarkerFaceColor',Colord(4,:));
% plot(real(k3)/k1 , imag(k3)/k1, 'd', 'markersize',4,...
%     'MarkerFaceColor',Colord(5,:));
% plot(real(Real_roots)/k1 , imag(Real_roots)/k1, 's', 'markersize',6,...
%     'MarkerFaceColor',Colord(5,:));
xlabel('$\Re\textrm{k}_{\rho}$','interpreter','latex')
ylabel('$\Im\textrm{k}_{\rho}$','interpreter','latex')
legend('Poles','Branch Point',...
    'Location','southwest','Orientation','horizontal');
if nu == 0
    title(['TE Pole Locations for thickness d = ', num2str(d2), 'm']);
else
    title(['TM Pole Locations for thickness d = ', num2str(d2), 'm']);
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

% Save tikz figure
% cleanfigure();
% if nu == 0
%     matlab2tikz('filename',sprintf('figures/TE_pole_loc_d_%d.tex',floor(d*1e7)),'showInfo', false);
% else
%     matlab2tikz('filename',sprintf('figures/TM_pole_loc_d_%d.tex',floor(d*1e7)),'showInfo', false);
% end



%% Plot Relative Pole Locations from the Branch Point
figure(2)
N = 2; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren');
Colord = get(gca, 'ColorOrder');
%
plot(real(root) , (imag(root)), 's', 'markersize',4,...
    'MarkerFaceColor',Colord(1,:));
hold on
plot(real(k4) , imag(k4) , 'd', 'markersize',4,...
    'MarkerFaceColor',Colord(2,:));
%
xlabel('$\textrm{Real Relative Distance from Branch Point}$','interpreter','latex')
ylabel('$\textrm{Imaginary Relative Distance from Branch Point}$','interpreter','latex')
legend('Poles','Branch Point',...
    'Location','northeast','Orientation','horizontal');
if nu == 0
    title(['Relative TE Pole Locations for thickness d = ', num2str(d2), 'm']);
else
    title(['Relative TM Pole Locations for thickness d = ', num2str(d2), 'm']);
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

% Save tikz figure
% cleanfigure();
% if nu == 0
%     matlab2tikz('filename',sprintf('figures/TE_pole_rel_loc_d_%d.tex',floor(d*1e7)),'showInfo', false)
% else
%     matlab2tikz('filename',sprintf('figures/TM_pole_rel_loc_d_%d.tex',floor(d*1e7)),'showInfo', false)
% end

%% Plot Pole Verification
figure(3)
N = 2; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% Obtain the colors of the plot to apply it on marker faces
Colord = get(gca, 'ColorOrder');
%
plot(real(D(root)), 's', 'markersize',4,...
    'MarkerFaceColor',Colord(1,:));
hold on
plot(imag(D(root)), 's', 'markersize',4,...
    'MarkerFaceColor',Colord(2,:));
%
xlabel('$\textrm{Real Relative Distance from Branch Point}$','interpreter','latex')
ylabel('$\textrm{Imaginary Relative Distance from Branch Point}$','interpreter','latex')
legend('Real Part','Imaginary Part',...
    'Location','southeast','Orientation','horizontal');
title('Evaluation of Denominator at pole locations');
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

% Save tikz figure
% cleanfigure();
% if nu == 0
%     matlab2tikz('filename',sprintf('figures/TE_pole_discrepancy_d_%d.tex',floor(d*1e7)),'showInfo', false)
% else
%     matlab2tikz('filename',sprintf('figures/TM_pole_discrepancy_d_%d.tex',floor(d*1e7)),'showInfo', false)
% end
% End
toc


