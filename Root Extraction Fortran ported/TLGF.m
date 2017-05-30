function D = TLGF(kp)
% ***********************************************************************
% 
%      Computes the Transmission Line Greens Function
% 
% ***********************************************************************

% Data from datafile of test_case3 from the reference paper
f = 1e12;
omega = 2*pi*f;


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

kz1 =  sqrt(k1 ^2 - kp .^2);

kz2 =  sqrt(k2 ^2 - kp .^2);

kz3 =  sqrt(k3 ^2 - kp .^2);

kz4 =  sqrt(k4 ^2 - kp .^2);

% Enforce kzn on the top sheet
for i = 1 : length(kz1)
    if imag(kz1(i)) >= 0
        kz1(i) = conj(kz1(i));
    end
    
%     if imag(kz2(i)) >= 0
%         kz2(i) = conj(kz2(i));
%     end
%     
%         if imag(kz3(i)) >= 0
%         kz3(i) = conj(kz3(i));
%     end
%     
    if imag(kz4(i)) >= 0
        kz4(i) = conj(kz4(i));
    end
end

% Define impedances
if nu == 0
    
    % TE Case
    Z1 =  omega./kz1;
    Z2 =  omega./kz2;
    Z3 =  omega./kz3;
    Z4 =  omega./kz4;
else
    % TM case
    Z1 =  kz1./(omega*ep1);
    Z2 =  kz2./(omega*ep2);
    Z3 =  kz3./(omega*ep3);
    Z4 =  kz4./(omega*ep4);
end


Gamma_12 =  (Z1 - Z2)./(Z1 + Z2);
Gamma_32 =  (Z3 - Z2)./(Z3 + Z2);
Gamma_43 =  (Z4 - Z3)./(Z4 + Z3);


% % Reflection Coefficients
Gamma_left =  (Gamma_12 + (-1)*exp(-2i*kz1*d1))...
    ./(1 + Gamma_12*(-1).*exp(-2i*kz1*d1)); % Left-looking
Gamma_right =  (Gamma_32 + Gamma_43.*exp(-2i*kz3*d3))...
    ./(1 + Gamma_32.*Gamma_43.*exp(-2i*kz3*d3)); % Right-looking

% Unknown A
A =  (Gamma_left .* exp(1i*kz2*2*z2))./(1 - Gamma_left.*Gamma_right.*exp(-2i * kz2 * d2)) ...
    .*( exp(-1i * kz2 * zp) + Gamma_right .* exp( -1i * kz2 * (2*z3 - zp)));

% Unknown B
B =  (Gamma_right .* exp(-1i*kz2*2*z3))./(1 - Gamma_left.*Gamma_right.*exp(-2i * kz2 * d2)) ...
    .*( exp(+1i * kz2 * zp) + Gamma_left .* exp( +1i * kz2 * (2*z2 - zp)));

% Denominator
D =  1 - Gamma_left.*Gamma_right.*exp(-2i * kz2 * dpp).*exp(-2i * kz2 * dp);

end