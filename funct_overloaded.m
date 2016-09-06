<<<<<<< HEAD
function y = funct_overloaded(varargin)
=======
function y = funct(varargin)
>>>>>>> master
% This function implements integrand of I_1(\rho)

% Courtesy of Mazin M Mustafa

switch nargin
    
    
    % Two arguments in the input. This is the normal case
    case 2
        c = varargin{1};
        d = varargin{2};
        x = c + d;
        % When singularity lies ahead
        if d >= 0
            y = ((1+x)^-0.5)*((1-x)^-0.5);
        else
            y = ((1+x)^-0.5)*((-d)^-0.5);
        end
        
        % When singularity lies behind
        %     if d <= 0
        %         y = ((1+x)^-0.5)*((1-x)^-0.5);
        %     end
        %     if d > 0
        %         y = (1/sqrt(-1))*((1+x)^-0.5)*((d)^-0.5);
        %     end
        % end
        % p = pi/2;
        % p = 53*pi/2;
        % y = besselj(0,p*x)*x*y;
        % % Examples from
        % y = x/(1+x^2)*besselj(0,x);
        % y = x^2*besselj(0,x);
        % y = 1/2*log(x^2 + 1)*besselj(1,x);
        y = (1-exp(-x))/(x*log(1 + sqrt(2)))*besselj(0,x);
        
        % Third argument will be used as sweep variable
    case 3
        c = varargin{1};
        d = varargin{2};
        z = varargin{3}; % Sweep variable
        x = c + d;
        % When singularity lies ahead
        if d >= 0
            y = ((1+x)^-0.5)*((1-x)^-0.5);
        else
            y = ((1+x)^-0.5)*((-d)^-0.5);
        end
        
        % When singularity lies behind
        %     if d <= 0
        %         y = ((1+x)^-0.5)*((1-x)^-0.5);
        %     end
        %     if d > 0
        %         y = (1/sqrt(-1))*((1+x)^-0.5)*((d)^-0.5);
        %     end
        % end
        
        % from reference paper [1] eq. 80
        p = 0;
        t = 1;
        y = exp(-x*z)*besselj(0,p*x)*besselj(3/2,x*t);
    otherwise
        error('Unexpected number of arguments');
end
end

% % Bessel Function reference
% @article{lucas1995evaluating,
%   title={Evaluating infinite integrals involving Bessel functions of arbitrary order},
%   author={Lucas, SK and Stone, HA},
%   journal={Journal of Computational and Applied Mathematics},
%   volume={64},
%   number={3},
%   pages={217--231},
%   year={1995},
%   publisher={Elsevier}
% }