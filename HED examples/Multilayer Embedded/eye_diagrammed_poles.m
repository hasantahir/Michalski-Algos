f = @(z) z.^50 + z.^12 - 5*sin(20*z).*cos(12*z) - 1;
lxlim = -20.3;
uxlim = 20.7;
kp = linspace(lxlim,uxlim,1e4);
f = chebfun(D,[lxlim uxlim],'splitting','on');
r = roots(f,'all','recursion');
plot(real(r), imag(r),'o')