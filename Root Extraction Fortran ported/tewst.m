clear;close all


% We declare the power of z as a global variable as ...
% we don't want to pass it on to the function as an argument. ...
% The Function itself is called through its handle
global pow 
%% Contour Box Definition
%% NEVELS MICHALSKI SPP PAPER
lambda = 633e-9; % Red light wavelength
eps_silver =  -18.295 - 1i*0.48085; % Johnson & Christy,1972 (refractiveindex.info) at 633 nm
% eps_silver = -265.06 - 1i*29.436;  % @ 2500 nm
load em_constants.mat % Contains varepsilon, mu and c
eps0 = epsilon_0;
mu0  = mu_0;
eps1 = 1 ;
c = 1/sqrt(mu0*eps0);
omega = 2*pi*c/lambda; % angular frequency
eta_0 = sqrt(mu0/eps0);
k_air = omega*sqrt(mu0*eps0); % propagation constant of air
k_silver = omega * sqrt(mu0*eps0*eps_silver); % propagation constant of silver


%% TLGF
f = 1e12;
omega = 2*pi*f;
lambda = 3e8/f;
num = 20; %Size of the arrays
tol = 1e-7;

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

% Contour Definition
% lxlim = 1*k2;
% uxlim = 1.00005*k2;
point = .1*k4 - .01i*k2;
step = 1.45*k4 + .02i*k2;

%% For dellnitz function with A B C T
% point = -1500 -1500i;
% step = 3000 + 3000i;

%% Gold Film example
% point = 1.00 - 1i*1.50;
% step = 1.05 + 1i*.20;

%% Polynomial


% point = -.20 -.20i;
% step = .40 + .40i;

npts = 8192;

%% Define the contour
%
% .... Lower side of rectangle
%

zstart = point;
zpoint_br = linspace(zstart, zstart + real(step), npts);
%
% .... Right hand vertical side of rectangle
%

zstart = point + real(step);
zpoint_ru = linspace(zstart, zstart + 1i*imag(step), npts);
%
% .... Top side of rectangle
%

zstart = point + step;
zpoint_ul = linspace(zstart, point + 1i*imag(step), npts);
%
% .... Left-hand vertical side of rectangle back to an including the
%      starting point
%

zstart = point + 1i*imag(step);
zpoint_lb = linspace(zstart, point, npts);


%% Concatenate all of them
% Last points are common so delete them

zpoint = [zpoint_br(1:end), zpoint_ru(1:end-1), zpoint_ul(1:end-1), zpoint_lb(1:end-1)];


%% Find roots based on the argument principal method
% Each 2*pi transition around the contour indicates presence of a zero
F = FZ(zpoint);
real_F = real(F);
imag_F = imag(F);
arg_F = atan2(imag_F, real_F);
%
% ---- Now try to count the transitions and bias the arguments
%
%      For example the function argument will jump from being
%      nearly pi to being nearly -pi.
%
%      We assign to each point a "bias" index. We don't want to
%      correct the arguments while we do this as it requires quite a
%      lot of book keeping.
%
%      At the end, we then, process all arguments and bias them
%      appropriately.
%
jbias = 1;
%
bias(1) = jbias;
%
for i = 2 : (npts*4-3)
    %
    %
    % ....... There are two cases
    %
    %            1. Argument goes from +ve value near to pi to
    %               a negative value near to pi.
    %
    %            2. Argument goes from -ve values near to pi to
    %               a positive value near to pi.
    %
    if ( arg_F(i-1) > pi/2  && (arg_F(i)   < -pi/2) )
        jbias = jbias + 1;
    elseif( arg_F(i-1) < -pi/2 && arg_F(i)   >  pi/2  )
        jbias = jbias - 1;
    end
    %
    bias(i) = jbias;
end
%
for i = 2 : (npts*4-3)
    arg_F(i) = arg_F(i) + real( (bias(i) - 1))*2*pi;
end

%
% ---- Finally, compute the number of roots
%
delta  = ( arg_F(npts*4-3) - arg_F(1) ) / (2*pi);
%
nroots = round(delta);   % Round off to nearest integer
s = zeros(1,2*nroots+1);
%

%% Adaptive Gauss-Quadrature Performed on the closed contour
for pow = 1 : 2*nroots + 1
    C = [zpoint(npts) zpoint(2*npts-1) zpoint(3*npts-2)];
    s(pow) = integral(@myfun, zpoint(1), zpoint(1), 'Waypoints', C,...
        'RelTol', 1e-12, 'AbsTol', 1e-12);    
end
s = s/(2*pi*1i);

