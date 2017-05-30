function D = TLGFd(kp)
% ***********************************************************************
%
%      Computes the Transmission Line Greens Function of a three
%      layer structure with source located at the center of the middle
%      layer
%
% ***********************************************************************


%% Calls polee_search_MIM
% Data from datafile of test_case3 from the reference paper
lambda = 1550e-9;
c = 3e8;
omega = 2*pi*c/lambda;
% -143.49668392243762-I*9.517339564114454

% Example Validations

% Material Properties
% Material Properties
ep1 = -143.499668392243762;
ep2 = 1;
ep3 = -143.499668392243762;

% ep1 = -143.49668392243762 +1i*9.517339564114454;
% ep2 = 1;
% ep3 = -143.49668392243762 +1i*9.517339564114454;
% ep1 = 2.0736;
% ep2 = -16 - 2i;
% ep3 = 2.0736;

% EM constants
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;

% Propagations Constants
k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);
k3 = omega*sqrt(mu0*ep0*ep3);

% Middle Layer thickness
d = lambda/4;
% d = 50e-9;

% Source Location
% zp = -d/2;



% % 
% kz1 =  sqrt(k1 ^2 - k1^2 - kp .^2);
% 
% kz2 =  sqrt(k2 ^2 - k1^2 - kp .^2);
% 
% kz3 =  sqrt(k3 ^2 - k1^2 - kp .^2);

% kz1 =  sqrt(k1 ^2 - kp .^2);
% 
% kz2 =  sqrt(k2 ^2 - kp .^2);
% 
% kz3 =  sqrt(k3 ^2 - kp .^2);

% kz1 =  sqrtbr(k1 ^2  - k1^2 - kp .^2, pi/2);
% 
% kz2 =  sqrt(k2 ^2  - k1^2 - kp .^2);
% 
% kz3 =  -sqrtbr(k3 ^2  - k1^2 - kp .^2, pi/2);

% TE/TM switch
nu = 1;

kz1 = (sqrtbr(k1 ^2 - k1^2 - kp .^2, -pi/2));
% 
kz2 = (sqrt(k2 ^2 - k1^2 - kp .^2));
% 
kz3 = (sqrtbr(k3 ^2 - k1^2 - kp .^2, pi/2));

% Normalize the propagation constants
% kair= omega*sqrt(mu0*ep0);
% 
% kz1 = kz1/kair;
% kz2 = kz2/kair;
% kz3 = kz3/kair;

% Enforce kzn on the top sheet
% fixidx = find(imag(kz1) < 0 & real(kz1) < 0);
% kz1(fixidx) = -kz1(fixidx);
% % 
% fixidx = find(imag(kz3) < 0 & real(kz3) < 0);
% kz3(fixidx) = -kz3(fixidx);
% if imag(kz1) >= 0
%     kz1 = conj(kz1);
% end
% 
% 
% if imag(kz3) <= 0
%     kz3 = conj(kz3);
% end




% Define impedances
if nu == 0
    
    % TE Case
    Z1 =  mu0*omega./kz1;
    Z2 =  mu0*omega./kz2;
    Z3 =  mu0*omega./kz3;
else
    % TM case
    Z1 =  kz1./(omega*ep1*ep0);
    Z2 =  kz2./(omega*ep2*ep0);
    Z3 =  kz3./(omega*ep3*ep0);
end

% Normalize Impedances to free-space
% Z0 = sqrt(mu0/ep0);
% Z1 =  Z1/Z0;
% Z2 =  Z2/Z0;
% Z3 =  Z3/Z0;



%% Reflection Coefficients
Gamma_left = (Z3 - Z2) ./ (Z3 + Z2); % Left-looking
Gamma_right =  (Z1 - Z2) ./ (Z1 + Z2); % Right-looking

% Denominator, Dispersion relation
D =  1 - Gamma_left.*Gamma_right.*exp(-2i * kz2 * d);

end