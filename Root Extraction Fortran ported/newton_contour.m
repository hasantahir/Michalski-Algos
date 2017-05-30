function [z,kount] = newton_contour(F,Fprime,z)
   % Newton.  [z,kount] = newton(F,Fprime,z).
   % Start Newton's method for function F and derivative Fprime
   % at a scalar complex point z.  Return converged value z and
   % the iteration kount. 
      sqrteps = sqrt(eps(2));
      kount = 0;
      w = Inf;
      while abs(z-w) > sqrteps
         kount = kount+1;
         w = z;
         z = z - F(z)./Fprime(z);
      end
   end