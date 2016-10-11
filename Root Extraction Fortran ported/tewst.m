clear;close all
f = 1e12;
omega = 2*pi*f;
lambda = 3e8/f;
num = 1e2; %Size of the arrays
global pow
% % % % % Example Validations
% % % % 
% % % % % Material Properties
% % % % ep1 = 1; % Air
% % % % ep2 = 9.7 ; % GaN/AlGaN layers combined
% % % % ep3 = 11; % Silicon base
% % % % 
% % % % % EM constants
% % % % mu0 = 4*pi*1e-7;
% % % % ep0 = 8.854e-12;
% % % % 
% % % % % Propagations Constants
% % % % k1 = omega*sqrt(mu0*ep0*ep1);
% % % % k2 = omega*sqrt(mu0*ep0*ep2);
% % % % k3 = omega*sqrt(mu0*ep0*ep3);
% % % % 
% % % % 
% % % % % Middle Layer thickness
% % % % d = .5*lambda;
% % % % 
% % % % % Source Location
% % % % zp = -d/2;
% % % % 
% % % % % Layer Heights
% % % % z0 = -d;
% % % % z1 = 0;

%% Example From Chinese paper
f = 5e9;
omega = 2*pi*f;
lambda = 3e8/f;
num = 5e2; %Size of the arrays

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

% point = -20.3 -20.7i;
% step = 40.6 + 41.4i;

point = -3*k4*(1 + .001i);
step = 6*k4*(1 + .001i);

% % % % point = 0 - k2*.01i;
% % % % step = 2*k3 + .02i*k2;
npts = 5012;
dnpts1 = npts - 1;
zstart = point;
zpoint_br = linspace(zstart, zstart + real(step), npts);
delta = real(step)/dnpts1;
zdelta_lower = delta;
% !
% !.... Right hand vertical side of rectangle
% !
% !     From
% !
% !       (point(1) + step(1) , point(2))
% !
% !     To
% !
% !       (point(1) + step(1), point(2) + step(2))
% !
zstart = point + real(step);
zpoint_ru = linspace(zstart, zstart + 1i*imag(step), npts);
delta = imag(step)/dnpts1;
zdelta_right = 1i*delta;
% !
% !.... Top side of rectangle
% !
% !     From
% !
% !       (point(1) + step(1) , point(2) + step(2))
% !
% !     To
% !
% !       (point(1) , point(2) + step(2))
% !
zstart = point + step;
zpoint_ul = linspace(zstart, point + 1i*imag(step), npts);
delta = -real(step)/dnpts1;
zdelta_upper = delta;
% !
% !.... Left-hand vertical side of rectangle back to an including the
% !     starting point
% !
% !     From
% !
% !       (point(1)  , point(2) + step(2))
% !
% !     To
% !
% !       (point(1) , point(2) )
% !
zstart = point + 1i*imag(step);
zpoint_lb = linspace(zstart, point, npts);
delta = -imag(step)/dnpts1;
zdelta_left = 1i*delta;

% Conctenate all of them
% Last points are common so delete them
zpoint = [zpoint_br(1:end), zpoint_ru(1:end-1), zpoint_ul(1:end-1), zpoint_lb(1:end-1)];
% !
% !---- Now compute the function and its argument all the way
% !     around the contour
% !
% !     Note that the arguments will normalized to the range
% !     [-pi, +pi] because of the behaviour of the ATAN2 function.
% !
F = FZ(zpoint);
real_F = real(F);
imag_F = imag(F);
arg_F = atan2(imag_F, real_F);
% !
% !---- Now try to count the transitions and bias the arguments
% !
% !     For example the function argument will jump from being
% !     nearly pi to being nearly -pi.
% !
% !     We assign to each point a "bias" index. We don't want to
% !     correct the arguments while we do this as it requires quite a
% !     lot of book keeping.
% !
% !     At the end, we then, process all arguments and bias them
% !     appropriately.
% !
jbias = 1;
% !
bias(1) = jbias;
% !
for i = 2 : (npts*4-3)
    % !
    % !
    % !....... There are two cases
    % !
    % !           1. Argument goes from +ve value near to pi to
    % !              a negative value near to pi.
    % !
    % !           2. Argument goes from -ve values near to pi to
    % !              a positive value near to pi.
    % !
    if ( arg_F(i-1) > pi/2  && (arg_F(i)   < -pi/2) )
        jbias = jbias + 1;
    elseif( arg_F(i-1) < -pi/2 && arg_F(i)   >  pi/2  )
        jbias = jbias - 1;
    end
    % !
    bias(i) = jbias;
end
% !
for i = 2 : (npts*4-3)
    arg_F(i) = arg_F(i) + real( (bias(i) - 1))*2*pi;
end

% !
% !---- Finally, compute the number of roots!
% !
delta  = ( arg_F(npts*4-3) - arg_F(1) ) / (2*pi);
% !
nroots = round(delta);   % Round off to nearest integer
s = zeros(1,2*nroots+1);
% !

%% Adaptive Gauss-Quadrature Performed on the closed contour
for pow = 1 : 2*nroots + 1
    %  !
    %  !....... Along lower side of rectangle
    %  !
    lower = quadgk(@myfun,zpoint(1),zpoint(npts),...
        'RelTol', 1e-7, 'AbsTol', 1e-7);
    %  !
    %  !....... Along right side of rectangle
    %  !
    right = quadgk(@myfun,zpoint(npts),zpoint(2*npts-1),...
        'RelTol', 1e-7, 'AbsTol', 1e-7);
    %  !
    %  !....... Along upper side of rectangle
    %  !
    top = quadgk(@myfun,zpoint(2*npts-1),zpoint(3*npts-2),...
        'RelTol', 1e-7, 'AbsTol', 1e-7);
    %  !
    %  !....... Along left side of rectangle
    %  !
    left = quadgk(@myfun,zpoint(3*npts-2),zpoint(4*npts-3),...
        'RelTol', 1e-7, 'AbsTol', 1e-7);
    s(pow) = (lower + right + top + left)/(2*pi*1i);
end
s = s/(2*pi*1i);

