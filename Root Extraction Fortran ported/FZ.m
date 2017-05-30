function f = FZ(z)
  % !***********************************************************************
  % !
  % !     Computes the function "f" at the point "z".
  % !
  % !***********************************************************************
  % !
  % !---- Tell the compiler to use the module which computes the dispersion
  % !     function.
  % !
  % use dispersion_function_module
  % !
  % !---- All variables in this unit have explicit type
  % !
%% Test Case 1
% f = exp(3*z) + 2*z.*cos(z) - 1;

%% Wilkinson Polynomial
% f =  1;
% for  ipoly = 1: 20
%     f = f .* (z - complex(ipoly));
% end

%% Modified Wilkinson Polynomial
f =  1;
for  ipoly = 1: 20
    f = f .* (z - complex(ipoly));
end
% Modified Wilkinson Polynomial with repeated roots
f = f .* ( z - 5).*(z - 6).^2  ...
.*(z - 7).^3;

%% Dellnitz paper
% f = 54.0 + z.*(44.0 + z.*(20.0 - z.*(3.0-z)));
%     df = 44.0 + 40*z - 9*z.^2 + 4*z.^3;
%     D = kp;
%     D = log - kp + 1.195281 + 0.536353*1i;
%   D = kp.^(-0.5*kp) - kp + 2.195093245 -2.750498700*1i;
% f = z.^50 + z.^12 - 5*sin(20*z).*cos(12*z) - 1;
% df = 50*z.^49 + 12*z.^11 - 100*cos(20*z).*cos(12*z)  +...
%     60*sin(20*z).*sin(12*z);
% D = z.^2 - 1;
%% Dispersion Function
% f = dispersion_function(z);
% f = dispersion_michalski(z);

%% Transmission Line Greens Function
%
% 1. PEC Backing
% f = TLGF(z);
%
% 2. Open Dielectrics
% f = TLGFd(z);
%% Example from Dellnitz paper
% A = -0.19435;
% B = 1000.41;
% C = 522463;
% T = 0.005;
% 
% f = z.^2 + A*z + B*exp(-T*z) + C;
% df = 2*z + A - B*T*exp(-T*z);
% f = 1e12;
% omega = 2*pi*f;
% lambda = 3e8/f;
% 
%% Simple sin function
% f = tan(z);
% f = z.^3-1;
end
