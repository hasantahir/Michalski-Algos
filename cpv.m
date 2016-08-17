% Evaluates the CPV integral of a user supplied function phi(x)
% over (a,b) at the point x0 in (a,b)
% by a n-point Double Epxponential (DE) scheme.

% The integal defintion is integral with respect to x over (a,b ) of  phi(x) /(x-x0).

% Double Exponential transformation between x & u is given by
% x= (b+a)/2 + (b-a)/2 tanh(pi/2 sinh(u))
% x0= (b+a)/2 + (b-a)/2 tanh(pi/2 sinh(u0))
% x, x0 lie in (a,b)
% u, u0 lie in (-inf,inf)
% Credits go to Euler, Mclaurin, C.Schwartz, Iri, Moriguti, Takasawa,
% Takahasi and Mori
function sum=cpv(a,b,x0,n)
format long;
% Given x0, invert to get u0;
arg=(2*x0-(b+a))/(b-a);u0=asinh((2/pi)*atanh(arg));
zmax=3.5;
h=2*zmax/n;
% Sampling points are chosen such that u0 lies between two successive
% quadrature points; Integral multiple of "h" is sliced off from (-zmax,zmax)
% and trapezoidal summation is done over over that interval, (umin,umax1)
[umax1,hu,umin,hl]=inter(zmax,u0,h);
z=umin-h;sum=0;
      while z <= umax1-h
      [x,dvdu]=de(a,b,z);
      arg=(phi(x))/(x-x0);
      der=dvdu*2*(x-a)*(b-x)/(b-a);
	   sum=sum+arg*der;
      z=z+h;
      end;
   sum=h*sum;
   %========================
   % Weight has to be (h/2) for z=umin & z=umax1
      z=umin;
      [x,dvdu]=de(a,b,z);
      arg=(phi(x))/(x-x0);
      der=dvdu*2*(x-a)*(b-x)/(b-a);
      sum=sum+der*arg*h/2;
      %==============================
      z=umax1;
      [x,dvdu]=de(a,b,z);
      arg=(phi(x))/(x-x0);
      der=dvdu*2*(x-a)*(b-x)/(b-a);
      sum=sum+der*arg*h/2;

% Add the contribution from the right fractional step by taking
% z as the midpoint of  -zmax & umin;
      z= (umin-zmax)/2;
      [x,dvdu]=de(a,b,z);
      arg=(phi(x))/(x-x0);
      der=dvdu*2*(x-a)*(b-x)/(b-a);
      sum=sum+hl*der*arg;
% Add the contribution from the right fractional step by taking
% z as the midpoint of  zmax & umax1;
      z=(zmax+umax1)/2;
      [x,dvdu]=de(a,b,z);
      arg=(phi(x))/(x-x0);
      der=dvdu*2*(x-a)*(b-x)/(b-a);
      sum=sum+hu*der*arg;

% cpv(a,b,x0,n) coding ENDS;