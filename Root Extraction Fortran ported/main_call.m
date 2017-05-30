% function main_call()
%  ---- Here we declare the namelist /INPUT/ that defines the search parameters.
%
%       It is useful to remember that
%
%
%          point  defines the complex value of the lower left hand co-ordinate
%                    of the rectangle
%
%          step   defines the length of the rectangle on the real and imaginary
%                    axes.
%
%          droot_converged  We declare that a root has been found whenever
%                           ABS(f(z)) < droot_converged
%
%          npts      Number of points in per side of rectangle in the fixed
%                    integration quadrature
%
%
clear all;clc;tic;close all
%% PEC Example
% 
% f = 1e12;
% omega = 2*pi*f;
% lambda = 3e8/f;
% % Material Properties
% ep1 = 12;
% ep2 = 4 ;
% ep3 = 2.1;
% ep4 = 1;
% 
% 
% % EM constants
% mu0 = 4*pi*1e-7;
% ep0 = 8.854e-12;
% 
% % Propagations Constants
% k1 = omega*sqrt(mu0*ep0*ep1);
% k2 = omega*sqrt(mu0*ep0*ep2);
% k3 = omega*sqrt(mu0*ep0*ep3);
% k4 = omega*sqrt(mu0*ep0*ep4);

%% Open Dielectric

% f = 1e12;
% omega = 2*pi*f;
% lambda = 3e8/f;
% % Material Properties
% ep1 = 1; % Air
% ep2 = 9.7 ; % GaN/AlGaN layers combined
% ep3 = 1; % Silicon base
% 
% 
% % EM constants
% mu0 = 4*pi*1e-7;
% ep0 = 8.854e-12;
% 
% % Propagations Constants
% k1 = omega*sqrt(mu0*ep0*ep1);
% k2 = omega*sqrt(mu0*ep0*ep2);
% k3 = omega*sqrt(mu0*ep0*ep3);

% point = -.20 -.20i;
% step = .40 + .40i;

%% For Dellnitz function with A B C T
% point = .95*k_air -.05*k_air*1i;
% step = .25*k_air + .12*k_air*1i;

%% Test Case 1
% point = -2.2 - 1i*3.5;
% step = 2.5 + 1i*8;

%% Gold-film example
% point = 1.9 - 1i*1.5;
% step = 0.4 + 1i*.20;

%% Wilkinson Polynomial
point = 4.5 - 1i*.5;
step = 2.0 + 1i*1.0;


%% TLGF
% point = 0.3*k1 - .002i*k1;
% step = 0.6*k1 + .004i*k1;

% point = .1*k2 - k2*.1i;
% step = 2.5*k2 + k2*.1i;

%% Simple functions
% point = -2 - 2i;
% step = 4 + 4i;
%----------------------------------------------------------------------
%% Code paramters
%----------------------------------------------------------------------

% Convergence Criterion
droot_converged = 1e-12;

% Number of points on each side of the contour
% Increase it to ensure capturing all Riemann sheets
npts = 8192*8;

% Currently not used as splitting algorithm not implemented
maxboxes = 500;

% Limit the number of roots per box
% If this is greater, split
max_roots_per_box = 5;

% Maxium roots to be found by the routine
maxroots = 500;



%----------------------------------------------------------------------
%% ROOT COUNTING and CONTOUR INTEGRALS
%----------------------------------------------------------------------
% Count the number of roots within the rectangle now by examination
% of the change in the argument and also compute the integrals "s",
% defined in the paper, need to construct the matrix eigenvalue
% problem.
%
% The contour integral is computed through MATLAB's integral function
% using way-points to define the path of integration

[nroots, s] = countz (point, step, npts, maxroots);


%----------------------------------------------------------------------
%% BUILD HANKEL MATRICES
%----------------------------------------------------------------------
%
% Construct H and H1 from the computed contour integrals s_n
%

for k = 1 : nroots
    for l = 1 : nroots
        H1(k,l) = s(k + l);
        H(k,l) = s(k + l - 1 );
    end
end

%----------------------------------------------------------------------
%% SOLVE EIGENVALUE PROBLEM
%----------------------------------------------------------------------
% The eigenvalue problem looks like
%
% H - \lambda \time H1 = 0
%
% The eigenvalues, \lambda are the initial roots of the system
%
zinitial_roots = eig(H1,H,'chol');

%----------------------------------------------------------------------
%% USE NEWTON / HALLEY'S METHOD TO REFINE THE ROOTS
%----------------------------------------------------------------------

qpl = point;
%
qpt = point + step;
%
% Check if we have a converged solution at this stage
% for each root. If not, apply Newton's / Halley's method
%
temp = [];

%----------------------------------------------------------------------
%% Check location of roots
%----------------------------------------------------------------------

% Before refining the roots, ensure they are inside the constructed contour
% If found outside, delete
true_roots = 0; % true number of roots
for i = 1 : nroots
    check = insidebox(qpt,qpl,zinitial_roots(i));
    if check
        true_roots = true_roots + 1;
        temp = vertcat(temp, zinitial_roots(i));
    end
end
zinitial_roots = temp;
% %
% % Check the quality of initial roots
% %
zinitial_func = FZ(zinitial_roots);
%

%----------------------------------------------------------------------
%% Refine roots
%----------------------------------------------------------------------

for k = 1 : true_roots
    %
    % If initially obtained roots are good enough
    %
    if(abs(zinitial_func(k)) < droot_converged)
        
        zfinal_roots(k) = zinitial_roots(k);
        zfinal_func(k)  = zinitial_func(k);
        %
        % Otherwise call Halley's method
        %
    else
        zfinal_roots(k) = halley3(qpl, qpt, zinitial_roots(k), droot_converged,true);
%         zfinal_roots(k) = newtzero(@FZ, zinitial_roots(k),50,1e-13);
% %         %
        % Compute the value of the function at the converged root.
        %
%         zfinal_func(k) = FZ(zfinal_roots(k));
    end
end

%----------------------------------------------------------------------
%% If we wish to find ALL roots
%----------------------------------------------------------------------

zfinal_roots = [];
for i = 1 : true_roots
    r = newtzero(@FZ, zinitial_roots(i));
    zfinal_roots = vertcat(zfinal_roots,r);
end
% Sort the array
zfinal_roots = sort(zfinal_roots);
% Clean up roots by weeding out too close values
if ~isempty(zfinal_roots)
    cnt = 1;  % Counter for while loop.
    
    while ~isempty(zfinal_roots)
        vct = abs(zfinal_roots - zfinal_roots(1)) < 1e-1; % Minimum spacing between roots.
        C = zfinal_roots(vct);  % C has roots grouped close together.
        [idx,idx] = min(abs(FZ(C)));  % Pick the best root per group.
        rt(cnt) = C(idx); %  Most root vectors are small.
        zfinal_roots(vct) = []; % Deplete the pool of roots.
        cnt = cnt + 1;  % Increment the counter.
    end
    zfinal_roots = sort(rt).';  % return a nice, sorted column vector
end
% clear zfinal_func zfinal_roots
zfinal_func = FZ(zfinal_roots)'
zfinal_roots = zfinal_roots'
toc
