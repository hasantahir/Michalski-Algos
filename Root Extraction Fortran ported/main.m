clear; clc
%  ---- Here we declare the namelist /INPUT/ that defines the search parameters.
%  
%       It is useful to remember that
%  
%  
%          point(1)  defines the REAL value of the lower left hand co-ordinate
%                    of the rectangle
%  
%          point(2)  defines the IMAGINARY value of the lower left hand co-ordinate
%                    of the rectangle
%  
%          step(1)   defines the length of the rectangle on the REAL axis.
%  
%          step(2)   defines the length of the  rectangle on the IMAGINARY axis
%  
%          droot_converged  We declare that a root has been found whenever
%                           ABS(f(z)) < droot_converged
%  
%          npts      Number of points in per side of rectangle in the fixed
%                    integration quadrature
%  
%  
% f = 1e12;
% omega = 2*pi*f;
% lambda = 3e8/f;
% % Material Properties
% ep1 = 1; % Air
% ep2 = 9.7 ; % GaN/AlGaN layers combined
% ep3 = 11; % Silicon base
% 
% % EM constants
% mu0 = 4*pi*1e-7;
% ep0 = 8.854e-12;
% 
% % Propagations Constants
% k1 = omega*sqrt(mu0*ep0*ep1);
% k2 = omega*sqrt(mu0*ep0*ep2);
% k3 = omega*sqrt(mu0*ep0*ep3);

%% Example From Chinese paper
f = 50e9;
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

point = k4 - k4*.1i;
step = 3.5*k4 + k4*.1i;
maxroots = 20;
droot_converged = 1e-9;

% Import from Namelist
% point  = -2.2 - 3.5*1i;
% step   = 5.0 + 8.0 *1i;

% point = -5 -5i;
% step = 10 + 10i;
% droot_converged = 1.0d-1;
zaitken = true;
zadaptive = false;
npts = 1024;
maxboxes = 500;
max_roots_per_box = 5;
% 
% point = -20.3 -20.7i;
% step = 40.6 + 41.4i;

% point = -k2 - k2*1i;
% step = 2*k2 + k2*2i;
%  
%  ----------------------------------------------------------------------
%  
%       BOX SPLITTING FOR THIS INPUT DOMAIN
%  
%  ----------------------------------------------------------------------
%  
%  ---- For multiple input sets, we'd bettwe deallocate any previously
%       allocated dynamic arrays for the box splitting.
%  
%  
%  ---- Initializes all boxes to have -1 roots. We use this as flag
%       later; remember that the splitting may produce some sub-boxes
%       with zero roots and that soem of the allocated space may not
%       even be used.
%  
% nroots = -1*ones(1, maxboxes);
%  
%  ---- Now do the work of splitting the contour
%  
% call  split_contour(point,step,&
%                     max_roots_per_box,points,steps,nrootsy,&
%                     maxboxes,'TRUE')
%  
% [points, step, nroots] = split_contour(point,step,max_roots_per_box,points,step,nroots, maxboxes,npts);
% [points, step, nroots, maxboxes,npts] =  split_contour(point,step,max_roots_per_box, npts);
% for ibox = 1 : maxboxes
%     if nrootsy <= 0
%         continue;
%     end
    %  
    %  ----------------------------------------------------------------------
    %  
    %       Step 1:  ROOT COUNTING and CONTOUR INTEGRALS
    %  
    %  ----------------------------------------------------------------------
    %  
    %  ---- Count the number of roots within the rectangle now by examination
    %       of the change in the argument and also compute the integrals "s",
    %       defined in the paper, need to construct the matrix eigenvalue
    %       problem.
    %  
    %       The integrals will becomputed wither using fixed point quadrature
    %       or by using an adaptive Gauus-Kronod scheme from QUADPACK.
    %  
    % call countz(points(1,ibox), steps(1,ibox), zadaptive, npts, maxroots, nroots, s)
    [nroots, s] = countz (point, step, npts, maxroots);
    %  
    %  ---- We'd better be careful that "nroots" does not
    %       exceed "maxroots".
    %  
    %       Also, we'd better check that actually have
    %       some roots to find.
    %  
    %  
    %  ---- Ok, so now we initialize the allocated space
    %  
    k = [1, nroots];
    l = [1, 8*nroots];
    zfinal_roots  =  zeros(k);
    zfinal_func  =  zeros(k);
%     radii         =  zeros(k);
%     indexv        =  zeros(k);
    vl1 = zeros(k(2)); % k by k array
    vr1 = zeros(k(2)); % k by k array
    work = zeros(l);
    rwork = zeros(l);
    %  
    %  ----------------------------------------------------------------------
    %  
    %       Step 2:  GENERALIZED EIGENVALUE PROBLEM
    %  
    %  ----------------------------------------------------------------------
    %  
    %  ---- Now build the matrices HMAT1 and HMAT as defined in the
    %       theory. The roots will be the eigenvalues, "lambda", of
    %       the system
    %  
    %           HMAT1 - lambda*HMAT
    %  
    for k = 1 : nroots
        for l = 1 : nroots
            H1(k,l) = s(k + l);
            H(k,l) = s(k + l - 1 );
        end
    end
    
    
    % Solve Eigenvalue problem
    %
    zinitial_roots = eig(H1,H, 'qz');
    %  
    %  ---- Solve the generalized eigenvalue problem using the LaPack routine
    zinitial_func = FZ(zinitial_roots);
    %  
    %  ----------------------------------------------------------------------
    %  
    %       Step 3:  ITERATE ROOTS USING HALLEY'S METHOD
    %  
    %  ----------------------------------------------------------------------
    %  
    %  ---- Prepare the co-ordinates of the rectangle in complex form
    %  
    qpl = point;
    %  
    qpt = point + step;
    %  
    %  ---- Ok, so let's see if we have a converged solution at this stage
    %       for each root. If not, apply Halley's method
    %  
    for k = 1 : nroots
        if(abs(zinitial_func(k)) < droot_converged)
            
            zfinal_roots(k) = zinitial_roots(k);
            zfinal_func(k)  = zinitial_func(k);
        else
            zfinal_roots(k) = halley3(qpl, qpt, zinitial_roots(k), droot_converged,zaitken);
            %  
            %  ......... Compute the value of the function at the converged root.
            %  
            zfinal_func(k) = FZ(zfinal_roots(k));
        end
    end
    %  
    %  ----------------------------------------------------------------------
    %  
    %       Step 4:  PREPARE FINAL RESULTS
    %  
    %  ----------------------------------------------------------------------
    %  
    %  ---- Compute an index such that the compted roots will be ordered
    %       in terms of increasing magnitudes (i.e radius from origin).
    %  
    radii = abs(zfinal_roots);
    indexv = indexx(nroots,radii);
% end
