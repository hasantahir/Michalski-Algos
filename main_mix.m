% This program and the functions used within performs integration using the
% Automatic tanh-sinh quadrature, also known as Double Exponential (DE)
% rule. The programs follow the Algorithm 2 described in [1].

% Author : Hasan Tahir Abbas
% Date: 08-09-2016

% [1]. Krzysztof A. Michalski & Juan R. Mosig (2016) Efficient computation of
%      Sommerfeld integral tails – methods and algorithms, Journal of Electromagnetic Waves and
%      Applications, 30:3, 281-317


% Define lower limit
% a = 1;
a = 0; % Integral I_2(5.13562, 0, 1)

% Define tolerance
tol =1e-15;

% Call the mixed quadrature routine
Sum = MixedQuad(a, tol)
