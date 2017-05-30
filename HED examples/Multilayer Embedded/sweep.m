function roots = sweep(f)

omega = 2*pi*f;
lambda = 3e8/f;
num = 15; %Size of the arrays
tol = 1e-8;
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

% % TE/TM switch
nu = 1;

kz1 = @(kp) sqrt(k1 ^2 - k4 ^2  - kp .^2);

kz2 = @(kp) sqrt(k2 ^2 - k4 ^2  - kp .^2);

kz3 = @(kp) sqrt(k3 ^2 - k4 ^2  - kp .^2);

kz4 = @(kp) sqrtbr(k4 ^2 - k4 ^2 - kp .^2);

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
D = @(kp) 1 - Gamma_left(kp).*Gamma_right(kp).*exp(-2i * kz2(kp) * dpp).*exp(-2i * kz2(kp) * dp);

lxlim = 1*k4;
uxlim = 3.5*k4;


%%%%%
%%
% Material Properties
% ep1 = 1; % Air
% ep2 = 9.7 ; % GaN/AlGaN layers combined
% ep3 = 11; % Silicon base
% % backed by a PEC layer
% 
% % EM constants
% mu0 = 4*pi*1e-7;
% ep0 = 8.854e-12;
% 
% % Propagations Constants
% k1 = omega*sqrt(mu0*ep0*ep1);
% k2 = omega*sqrt(mu0*ep0*ep2);
% k3 = omega*sqrt(mu0*ep0*ep3);
% 
% 
% % Middle Layer thickness
% d = .5*lambda;
% 
% % Source Location
% zp = -d/2;
% 
% % Layer Heights
% z0 = -d;
% z1 = 0;
% z_pec = -5*d;
% 
% % TE/TM switch
% nu = 1;
% 
% kz1 = @(kp) sqrt(k1 ^2 - kp .^2);
% 
% kz2 = @(kp) sqrt(k2 ^2 - kp .^2);
% 
% kz3 = @(kp) sqrt(k3 ^2 - kp .^2);
% 
% % Define impedances
% if nu == 0
%     
%     % TE Case
%     Z1 = @(kp) omega./kz1(kp);
%     Z2 = @(kp) omega./kz2(kp);
%     Z3 = @(kp) omega./kz3(kp);
% else
%     % TM case
%     Z1 = @(kp) kz1(kp)./(omega*ep1);
%     Z2 = @(kp) kz2(kp)./(omega*ep2);
%     Z3 = @(kp) kz3(kp)./(omega*ep3);
% end
% % Gammas
% Gamma_32 = @(kp)(Z3(kp) - Z2(kp))./ (Z3(kp) + Z2(kp));
% Gamma_43 =  -1;
% 
% % % Left Looking Input Impedance
% % Z_in_left = @(kp)  Z2(kp) .* (Z3(kp) + Z2(kp)*1i.*tan(kz2(kp)*d/2))...
% %     ./(Z2(kp) + Z3(kp)*1i.*tan(kz2(kp)*d/2));
% % 
% % % Right Looking Input Impedance
% % Z_in_right = @(kp) Z2(kp) .* (Z1(kp) + Z2(kp)*1i.*tan(kz2(kp)*d/2))...
% %     ./(Z2(kp) + Z1(kp)*1i.*tan(kz2(kp)*d/2));
% % 
% % % % Reflection Coefficients
% % Gamma_left = @(kp) (Z_in_left(kp) - Z2(kp)) ./ (Z_in_left(kp) + Z2(kp)); % Left-looking
% % Gamma_right = @(kp)  (Z_in_right(kp) - Z2(kp)) ./ (Z_in_right(kp) + Z2(kp)); % Right-looking
% 
% % % Reflection Coefficients
% % Gamma_left = @(kp)(Z3(kp) - Z2(kp)) ./ (Z3(kp) + Z2(kp)); % Left-looking
% Gamma_left = @(kp)(Gamma_32(kp) + (-1)*exp(-1i*kz3(kp)*4*d))...
%     ./(1 + Gamma_32(kp).*(-1).*exp(-1i*kz3(kp)*4*d)); % Left-looking
% Gamma_right = @(kp) (Z1(kp) - Z2(kp)) ./ (Z1(kp) + Z2(kp)); % Right-looking
% 
% % Unknown A
% A = @(kp) (Gamma_left(kp) .* exp(1i*kz2(kp)*2*z0))./(1 - Gamma_left(kp).*Gamma_right(kp).*exp(-2i * kz2(kp) * d)) ...
%     .*( exp(-1i * kz2(kp) * zp) + Gamma_right(kp) .* exp( -1i * kz2(kp) * (2*z1 - zp)));
% 
% % Unknown B
% B = @(kp) (Gamma_right(kp) .* exp(-1i*kz2(kp)*2*z1))./(1 - Gamma_left(kp).*Gamma_right(kp).*exp(-2i * kz2(kp) * d)) ...
%     .*( exp(+1i * kz2(kp) * zp) + Gamma_left(kp) .* exp( +1i * kz2(kp) * (2*z0 - zp)));
% 
% % Denominator
% D = @(kp) 1 - Gamma_left(kp).*Gamma_right(kp).*exp(-2i * kz2(kp) * d);
% 
% lxlim = k1;
% uxlim = k3;
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
%         vct = abs(root - root(1)) < 5e-6; % Minimum spacing between roots.
%         C = root(vct);  % C has roots grouped close together.
%         [idx,idx] = min(abs(D(C)));  % Pick the best root per group.
%         rt(cnt) = C(idx); %  Most root vectors are small.
%         root(vct) = []; % Deplete the pool of roots.
%         cnt = cnt + 1;  % Increment the counter.
%     end
%     root = sort(rt).';  % return a nice, sorted column vector
% end
roots = root((real(root))>k4);% & abs(imag(root))<tol);
roots = roots/k4;
if ~isempty(roots)
    cnt = 1;  % Counter for while loop.
    
    while ~isempty(roots)
        vct = abs(roots - roots(1)) < 1e-2; % Minimum spacing between roots.
        C = roots(vct);  % C has roots grouped close together.
        [idx,idx] = min(abs(D(C)));  % Pick the best root per group.
        rt(cnt) = C(idx); %  Most root vectors are small.
        roots(vct) = []; % Deplete the pool of roots.
        cnt = cnt + 1;  % Increment the counter.
    end
    roots = sort(rt).';  % return a nice, sorted column vector
end
end