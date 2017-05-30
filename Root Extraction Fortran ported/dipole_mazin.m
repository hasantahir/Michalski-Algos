close all; clear ; clc;
fig = 1; j = sqrt(-1);
N = 1e2; % Maximum 1e2 NEVER exceed!
f = 1; w = 2*pi*f; c = 3e8; lamda = c/f; k = w/c;
L = input('Length of dipole as a fraction of lamda: '); L = L*lamda;
R = 2*lamda; x = linspace(-R,R,N); z = linspace(-R,R,N);
[xx,zz] = meshgrid(x,z);
r = sqrt(xx.^2 + zz.^2); theta = pi/2 - atan(zz./xx);
Nstep = 100; T = 2;
phase = linspace(0,2*pi*T,Nstep);
E = zeros(N,N,Nstep);
for n = 1:length(phase)   	
    E_theta = exp(j*phase(n))*exp(-j*k*r).*( (cos(0.5*k*L*cos(theta)) - cos(0.5*k*L))./(sin(theta)) );
    E(:,:,n) = real(E_theta);
end
figure(fig)
fig = fig + 1;
hold on
xx = xx./lamda; zz = zz./lamda; R = R/lamda;
for n = 1:Nstep 
    pcolor(xx,zz,E(:,:,n));
    shading interp
    colormap(brewermap([],'RdYlBu')) %RdYlGn 'RdYlBu' Spectral
    axis equal
    xlabel('x (\lambda)')
    ylabel('z (\lambda)')
    str = sprintf('E-field at \\phi = 0^{\\circ}, L = %0.2f \\lambda',L/lamda);
    title(str)
    axis([-R R -R R])
    pause(0.05);
end
hold off
N = 1e4;
theta = linspace(-pi,pi,N);
P = ( (cos(0.5*k*L*cos(theta)) - cos(0.5*k*L))./(sin(theta)) ).^2; P = P./max(max(P));
figure(fig)
fig = fig + 1;
plot(theta*180/pi,P,'k')
xlabel('\theta^{\circ}')
ylabel('U')
str = sprintf('Radiation pattern at \\phi = 0^{\\circ}, L = %0.2f \\lambda',L/lamda);
title(str)
xlim([-180 180])
figure(fig)
fig = fig + 1;
polar(theta,P,'k')
str = sprintf('Radiation pattern at \\phi = 0^{\\circ}, L = %0.2f \\lambda',L/lamda);
title(str)
view (-270,-90)