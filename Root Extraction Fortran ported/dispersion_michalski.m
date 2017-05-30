function out = dispersion_michalski(gamma)
% ***********************************************************************
% 
%      Computes the dispersion function defined in the  Nevels Michalski SPP
%      paper.
% 
%      Note that you must call the entry point
% 
%           set_lambda_for_dispersion_function()
% 
%      before this function is ever called. The entry point preloads
%      all necessary data whihc we then use at each and every point,
%      gamma, in the complex plane at which we compute the dispersion
%      function.
% 
% ***********************************************************************

% Data from datafile of test_case3 from the reference paper
lambda = 631.42; % in nanometers


epsilon1 = 1;
% epsilonm = -0.291844 +1i*6.511223;
% epsilonm = -10.5188  +1i*1.14385;
epsilonm = -18.295 - 1i*0.48085;

k1  = 2*pi/lambda;
k2  = 2*pi/lambda*sqrt(epsilonm);


% Arguments of the upcoming trig functions
gammasqd = gamma.^2;
% 

kz1 = sqrt(k1^2 - gammasqd);
kz2 = sqrt(k2^2 - gammasqd);

for i = 1 : length(kz1)
    if imag(kz1(i)) <= 0
        kz1(i) = conj(kz1(i));
    end
    
    if imag(kz2(i)) >= 0
        kz2(i) = conj(kz2(i));
    end
    
end
%
out = kz2/epsilonm + kz1/epsilon1;
kxp = k1*sqrt(epsilon1*epsilonm./(epsilon1 + epsilonm));
save kk.mat kxp
end