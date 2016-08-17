clear all; close all
nu = 1;
p = 1e1;
f = 10e9;
omega = 2*pi*f;
ep1 = 1;
ep2 = 10 - 1i*18;
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;
k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);
kz1 = @(x) sqrt(k1^2 - (x).^2);
kz2 = @(x) sqrt(k2^2 - (x).^2);
gamma_1h = @(x) -(kz2(x) - kz1(x))./(kz2(x) +kz1(x));
gamma_1e = @(x) (kz2(x)/ep2 - kz1(x))./(kz2(x)/ep2 +kz1(x));
G_1 = @(x) k1./(1i*kz1(x)).*(gamma_1h(x));
G_2 = @(x) k1./((x)).*(gamma_1e(x) - gamma_1h(x));
f = @(x) G_2(x).*besselj(1, p*(x)).*(x);
f1 = @(x) G_1(x).*besselj(0, p*(x)).*(x);
xx = linspace(1,20,21);
figure(1)
y = f(xx);
plot(xx,real(y))
hold on
plot(xx, imag(y))
hold off
p = linspace(1e-3/k1,1e1/k1,250);
figure(2)
y1 = f1(xx);
plot(xx,real(y1))
hold on
plot(xx, imag(y1))
hold off
p = linspace(1e-3/k1,1e1/k1,250);
a = 2*k1;
for i = 1 : length(p)
    

    f = @(x) G_2(x).*besselj(1, p(i)*(x)).*(x);
    f1 = @(x) G_1(x).*besselj(0, p(i)*(x)).*(x);

%     yy(i) = integral(f,0,.5e7,'RelTol',1e-12,'AbsTol',1e-12);
    yy1(i) = integral(f1,a,1e10,'RelTol',1e-12,'AbsTol',1e-12);
    warning('off');

end
% figure(3)
% loglog(k1*p,abs(yy)/k1)
figure(4)
loglog(k1*p,abs(yy1)/k1)