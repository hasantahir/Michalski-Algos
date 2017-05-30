function D = TLHEMT(kp)
% ***********************************************************************
%
%      Computes the Transmission Line Greens Function of a three
%      layer structure with source located at the center of the middle
%      layer
%
% ***********************************************************************


%% Calls polee_search_MIM
% Data from datafile of test_case3 from the reference paper
f = 1e12;
omega = 2*pi*f;
lambda = 3e8/f;
num = 80; %Size of the arrays

% Example Validations

% Material Properties
% ep1 = 1; % Air
% ep2 = 9.7 ; % GaN/AlGaN layers combined
% ep3 = 1; % Silicon base
% % backed by a PEC layer



% Chinese Homotopy Method
ep1 = 1; % Air
ep2 = 7.34 ; % GaN/AlGaN layers combined
ep3 = 11; % Silicon base
% EM constants
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;

% Propagations Constants
k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);
k3 = omega*sqrt(mu0*ep0*ep3);

% 
% if 
% kz1 =  sqrt(k1 ^2  - k1^2 - kp .^2);
% 
% kz2 =  sqrt(k2 ^2  - k1^2 - kp .^2);
% 
% kz3 =  sqrt(k3 ^2  - k1^2 - kp .^2);

% kz1 =  sqrtbr(k1 ^2  - k1^2 - kp .^2, pi/2);
% 
% kz2 =  sqrt(k2 ^2  - k1^2 - kp .^2);
% 
% kz3 =  -sqrtbr(k3 ^2  - k1^2 - kp .^2, pi/2);

% TE/TM switch
nu = 1;

kz1 = (1*sqrtbr(k1 ^2 - kp .^2, pi/2));

kz2 = (+sqrt(k2 ^2 - kp .^2));

kz3 = (-1*sqrtbr(k3 ^2 - kp .^2, -pi/2));

% Normalize the propagation constants
% kair= omega*sqrt(mu0*ep0);
% 
% kz1 = kz1/kair;
% kz2 = kz2/kair;
% kz3 = kz3/kair;

% Enforce kzn on the top sheet
% fixidx = find(imag(kz1) > 0 & real(kz1) > 0);
% kz1(fixidx) = conj(kz1(fixidx));
% % 
% fixidx = find(imag(kz3) < 0 & real(kz3) > 0);
% kz3(fixidx) = conj(kz3(fixidx));
% if imag(kz1) >= 0
%     kz1 = conj(kz1);
% end
% 
% 
% if real(kp) > k1 & real(kp) < k3
%     if imag(kz3) <= 0
%         kz3 = conj(kz3);
%     end
%     if imag(kz1) >= 0
%         kz1 = conj(kz1);
%     end
%     if real(kz3) > 0
%         kz3 = -conj(kz3);
%     end
%     if real(kz1) < 0
%         kz1 = -conj(kz1);
%     end
%     
% end
% 
% if real(kp) < k1 
%     if imag(kz3) <= 0
%         kz3 = conj(kz3);
%     end
%     if imag(kz1) <= 0
%         kz1 = conj(kz1);
%     end
%     if real(kz3) > 0
%         kz3 = -conj(kz3);
%     end
%     if real(kz1) > 0
%         kz1 = -conj(kz1);
%     end
%     
% end
% 
% if real(kp) > k3 
%     if imag(kz3) >= 0
%         kz3 = conj(kz3);
%     end
%     if imag(kz1) >= 0
%         kz1 = conj(kz1);
%     end
%     if real(kz3) < 0
%         kz3 = -conj(kz3);
%     end
%     if real(kz1) < 0
%         kz1 = -conj(kz1);
%     end
%     
% end




% Define impedances
if nu == 0
    
    % TE Case
    Z1 =  omega./kz1;
    Z2 =  omega./kz2;
    Z3 =  omega./kz3;
else
    % TM case
%     Z1 =  kz1./(omega*ep1*ep0);
%     Z2 =  kz2./(omega*ep2*ep0);
%     Z3 =  kz3./(omega*ep3*ep0);
    
    Z1 =  kz1./(omega*ep1);
    Z2 =  kz2./(omega*ep2);
    Z3 =  kz3./(omega*ep3);
end

% Normalize Impedances to free-space
% Z0 = sqrt(mu0/ep0);
% Z1 =  Z1/Z0;
% Z2 =  Z2/Z0;
% Z3 =  Z3/Z0;

% Middle Layer thickness
d = 200.2e-6;


%% Reflection Coefficients
Gamma_left = (Z3 - Z2) ./ (Z3 + Z2); % Left-looking
Gamma_right =  (Z1 - Z2) ./ (Z1 + Z2); % Right-looking

% Denominator, Dispersion relation
D =  1 - Gamma_left.*Gamma_right.*exp(-2i * kz2 * d);

end