% This program and the functions used within performs integration using the
% Automatic tanh-sinh quadrature, also known as Double Exponential (DE)
% rule. The programs follow the Algorithm described in [1].

% Author : Hasan Tahir Abbas
% Date: 08-09-2016

% [1]. Krzysztof A. Michalski & Juan R. Mosig (2016) Efficient computation of
%      Sommerfeld integral tails – methods and algorithms, Journal of Electromagnetic Waves and
%      Applications, 30:3, 281-317


% Define limits
a = 0;
b =1;

% Define tolerance
tol =1e-15;

% Call the tanh-sinh quadrature routine
Sum = TanhSinhQuad(a, b, tol)
