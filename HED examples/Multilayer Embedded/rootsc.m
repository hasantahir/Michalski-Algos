function w = rootsc(n, x, display)
%ROOTSC - Computes complex roots
% 
% w = ROOTSC(n,x,1) extracts the n n-th roots of the complex number x
%                   and (optionally) plots the roots themselves
%
% INPUTS:
% n: the root degree
% x: a complex number
% display: if display == 1 a final plot is provided with displays x and w
%
% OUTPUTS:
% w: a complex vector of length n, containing the n n-th roots of x
%
%
% If x = 1 and n = 2 than w = [-1, 1]
% If x = 1 and n = 3 than w = [-0.5 + 0.866i, -0.5 - 0.866i, 1]
% If x = 1 and n = 4 than w = [i, -1, -i, 1]
% ...
% 
% To check the attained results, verify that abs(w.^n-x*ones(1,n)) 
% is sufficiently small
%
% Example
% 
%       w = rootsc(5,(0.25+0.25*i),1); 
% 
%
% Please note that the same output w can be computed using 
% the built-in Matlab function ROOTS as follows:  
%
%       w = roots([1, zeros(1,n-1), -x]);
%
% However, ROOTSC can also plot the founded roots.
% Furthermore, it takes less computational time than ROOTS, 
% as the following example shows:
%
% Example
%
% k = 100;  % number of random complex numbers to use for test
% n = 20;   % the root degree
% 
% re = randn(k,1);
% im = randn(k,1);
% x = re + sqrt(-1).* im;
% 
% et_roots  = 0; % ellapsed time using ROOTS
% et_rootsc = 0; % ellapsed time using ROOTSC
% 
% for j=1:k
%     tic; w_roots  = roots([1, zeros(1,n-1), -x(j)]); et_roots = et_roots + toc;
%     tic; w_rootsc = rootsc(n,x(j),0); et_rootsc = et_rootsc + toc;
% end
% disp(sprintf('Ellapsed time using Matlab function ROOTS:  %g', et_roots));
% disp(sprintf('Ellapsed time using ROOTSC:                 %g', et_rootsc));
%
%
% 
% See also ROOTS, SQRT, SQRTM, REALSQRT, HYPOT.

% ROOTSC, $Version: 1.1, April 2006
% Author: Marco Cococcioni
% Please report any bugs to m.cococcioni <at> gmail.com
%
% History:
% Version 1.1: obsolete function PHASE replaced with function ANGLE,
%              inputs check added, second example added    
% Version 1.0: base version with no inputs check. This version uses the
%              obsolete function PHASE
% Acknowledgments: I would like to thank Mr. John D'Errico you its useful
%                  suggestions and remarks.


if nargin < 1,
    disp('Example of usage of ROOTSC:');
    disp(' ');
    disp('The call:')
    disp(' ');
    disp('ans = rootsc(5,(0.25+0.25*i),1)');
    disp(' ');
    disp('will compute the 5 complex roots of the complex number 0.25+0.25i');
    disp('The founded roots will be plotted in a polar diagram.');
    disp(' ');
    disp('Press a key to run rootsc with default values:');
    disp(' ');
    pause

    n = 5;
    x = 0.25+0.25*(sqrt(-1));
    display = 1;
elseif nargin < 2,
    x = 1; % computes the n roots of the unity
    display = 0;
elseif nargin < 3,
    display = 0;
end

if (not(floor(n)==n) || n <= 0),
    error('First argument n must be a positive integer number');
end


re_x = real(x);
im_x = imag(x);
angle_x = atan2(im_x, re_x);
abs_x   = sqrt(re_x*re_x + im_x*im_x);
ro = abs_x^(1/n);
due_pi = 2*pi;
delta_phi = due_pi/n;
phi_0 = angle_x/n;
phi = phi_0 + (delta_phi:delta_phi:due_pi);

w = ro*(cos(phi) + i* sin(phi));

if display,
    figure;
    if abs_x > ro,
        polar([0, angle_x], [0, abs_x], '-b.');
        set(gca,'nextplot','add');
        polar(angle(w([1:end,1])), abs(w([1:end,1])), ':m.');
    else
        polar(angle(w([1:end,1])), abs(w([1:end,1])), ':m.');
        set(gca,'nextplot','add');
        polar([0, angle_x], [0, abs_x], '-b.');
    end
    xlabel(sprintf('The %d complex roots of %g+%gi', n, real(x), imag(x)));
end