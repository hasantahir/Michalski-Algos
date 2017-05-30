function df = dFZ(z)
  % ***********************************************************************
  % 
  %      Computes the derivative of the function "f"
  %      at the point "z".
  % 
  % ***********************************************************************


sqrteps = (eps)^(1/3); % Used to make h.  
h = sqrteps*z; % From numerical recipes, make h = h(xr)
df = (FZ(z+h)-FZ(z-h))./(2*h); % First Derivative
% df = (-1/12*FZ(z+2*h)+2/3*FZ(z+h)-2/3*FZ(z-h) + 1/12*FZ(z-2*h))./(h); % First Derivative
% d2f = (FZ(z+h)+FZ(z-h) - 2*FZ(z))./(h.^2); % Second Derivative
end
