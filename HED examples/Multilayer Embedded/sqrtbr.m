function fval = sqrtbr(z, shiftbr)

% shifts the branch cut of square root function

% default branch shift cut
if nargin == 1
  shiftbr = -pi/2; % to the negative imaginary axis
end

fval = sqrt(z.*exp(-1i.*shiftbr)).*exp(-1i.*shiftbr./2);
