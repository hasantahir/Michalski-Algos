function [nroots, s] = countz (point, step, npts, maxroots)
  % SUBROUTINE countz(point, step, zadaptive, npts, maxroots, nroots, s)

  % !***********************************************************************
  % !
  % !     Computes the argument of the function f(z) around the contour
  % !     thereby determining the number of f(z) = 0 lying roots within.
  % !     Also computes the integrals "s"
  % !
  % !
  % !     Input data:
  % !        point(1) defines the REAL value of the lower left hand
  % !                 co-ordinate of the rectangular contour
  % !        point(2) defines the IMAGINARY value of the lower left hand
  % !                 co-ordinate of the rectangular contour
  % !        step(1)  defines the length of the rectangle on the REAL axis.
  % !        step(2)  defines the length of the rectangle on the IMAGINARY
  % !                  axis
  % !       zadaptive when .true. the QUADPACK implementation of an
  % !                 adaptive Gauss-Kronod scheme is used. Otherwise,
  % !                 fixed point trapezoidal is used.
  % !            npts the number of points to be used in the fixed point
  % !                 integration scheme - ignored if zadaptive=.true.
  % !        maxroots the dimension of s is s(0:2*maxroots)
  % !
  % !     Output data:
  % !          nroots the number of roots of f(z) lying within the contour
  % !               s the contour integrals of (powers of z)/f(z).
  % !
  % !     The routine terminates the program if the computed value for
  % !     nroots exceeds "maxroots".
  % !
  % !***********************************************************************
  global F % function value
  global pow
  % !
  % !---- Initialize the number of roots found to be zero, of course!
  % !
  % !     And of course the integrals returned.
  % !

  s = zeros(1,2*maxroots+1);

  % !
  % !---- Build the vector of points at which we shall evaluate the
  % !     function and therefore its argument.
  % !
  % !     The following figure displays the co-ords of each part
  % !     of the rectangle. We put "npts" points along each side.
  % !     Of course, the corner points count double - hence some
  % !     awkward indexing of the arrays.
  % !
  % !                                              Complex x plane
  % !
  % !
  % !   (point(1), point(2) + step(2))      (point(1) + step(1),
  % !                                             point(2) + step(2))
  % !       +----------------------------------------------------+
  % !       |                                                    |
  % !       |                                                    |
  % !       |                                                    |
  % !       |                                                    |
  % !       |                                                    |
  % !       |                                                    |
  % !       |                                                    |
  % !       |                                                    |
  % !       |                                                    |
  % !       +----------------------------------------------------+
  % !
  % !   (point(1),point(2))                 (point(1) + step(1), point(2))
  % !
  % !
  % !
  % !     We do this along each of the four sides of the rectangle
  % !     one at a time. We start6 the the lower left hand side and
  % !     then work anti-clockwise around all four sides.
  % !
  % !
  % !.... Lower side of rectangle
  % !
  % !        From
  % !
  % !             (point(1), point(2))
  % !
  % !        To
  % !
  % !             (point(1) + step(1), point(2))
  % !
  % Assuming point is complex number
  zstart = point;
  zpoint_br = linspace(zstart, zstart + real(step), npts);
  % !
  % !.... Right hand vertical side of rectangle
  % !
  % !     From
  % !
  % !       (point(1) + step(1) , point(2))
  % !
  % !     To
  % !
  % !       (point(1) + step(1), point(2) + step(2))
  % !
  zstart = point + real(step);
  zpoint_ru = linspace(zstart, zstart + 1i*imag(step), npts);
  % !
  % !.... Top side of rectangle
  % !
  % !     From
  % !
  % !       (point(1) + step(1) , point(2) + step(2))
  % !
  % !     To
  % !
  % !       (point(1) , point(2) + step(2))
  % !
  zstart = point + step;
  zpoint_ul = linspace(zstart, point + 1i*imag(step), npts);
  % !
  % !.... Left-hand vertical side of rectangle back to an including the
  % !     starting point
  % !
  % !     From
  % !
  % !       (point(1)  , point(2) + step(2))
  % !
  % !     To
  % !
  % !       (point(1) , point(2) )
  % !
  zstart = point + 1i*imag(step);
  zpoint_lb = linspace(zstart, point, npts);

  % Conctenate all of them
  % Last points are common so delete them
  zpoint = [zpoint_br(1:end), zpoint_ru(1:end-1), zpoint_ul(1:end-1), zpoint_lb(1:end-1)];
  % !
  % !---- Now compute the function and its argument all the way
  % !     around the contour
  % !
  % !     Note that the arguments will normalized to the range
  % !     [-pi, +pi] because of the behaviour of the ATAN2 function.
  % !
  F = FZ(zpoint);
  real_F = real(F);
  imag_F = imag(F);
  arg_F = atan2(imag_F, real_F);
  % !
  % !---- Now try to count the transitions and bias the arguments
  % !
  % !     For example the function argument will jump from being
  % !     nearly pi to being nearly -pi.
  % !
  % !     We assign to each point a "bias" index. We don't want to
  % !     correct the arguments while we do this as it requires quite a
  % !     lot of book keeping.
  % !
  % !     At the end, we then, process all arguments and bias them
  % !     appropriately.
  % !
  jbias = 1;
  % !
  bias(1) = jbias;
  % !
  for i = 2 : (npts*4-3)
    % !
    % !
    % !....... There are two cases
    % !
    % !           1. Argument goes from +ve value near to pi to
    % !              a negative value near to pi.
    % !
    % !           2. Argument goes from -ve values near to pi to
    % !              a positive value near to pi.
    % !
    if ( arg_F(i-1) > pi/2  && (arg_F(i)   < -pi/2) )
      jbias = jbias + 1;
    elseif( arg_F(i-1) < -pi/2 && arg_F(i)   >  pi/2  )
      jbias = jbias - 1;
    end
    % !
    bias(i) = jbias;
  end
  % !
  for i = 2 : (npts*4-3)
    arg_F(i) = arg_F(i) + real( (bias(i) - 1))*2*pi;
  end

  % !
  % !---- Finally, compute the number of roots!
  % !
  delta  = ( arg_F(npts*4-3) - arg_F(1) ) / (2*pi);
  % !
  nroots = round(delta);   % Round off to nearest integer
  % !
  s = zeros(1,2*nroots+1);
  % !
  % !---- Check "nroots" for validity
  % !
  % !     Exit this routine
  % !
  % !         If 0
  % !
  % !         If nroots > maxroots (the dimensions are insufficient)
  % !
  % if(nroots .le. 0)        goto 800
  % !
  % if(nroots .gt. maxroots) goto 800

  % !
  % !=======================================================================
  % !
  % !     C O M P U T E   "S"   I N T E G R A L S
  % !
  % !=======================================================================
  % !
  % !---- We now perform the integrals "s" defined in the paper
  % !
  % !     We have two options, adaptive or fixed rule. These
  % !     are implemented as two in-line subroutines to save
  % !     even more call overhead. Both in-line routines return
  % !     to line 700 on completion.
  % !
  % !     The fixed point scheme can re-use the computed function
  % !     values used above. The adaptibe does not need them, so we
  % !     can optimally free the memory at this point!
  % !
  % !
  %% Adaptive Gauss-Quadrature Performed on the closed contour
  for pow = 1 : 2*nroots + 1
    %  !
    %  !....... Along lower side of rectangle
    %  !
    lower = quadgk(@myfun,zpoint(1),zpoint(npts),...
    'RelTol', 1e-8, 'AbsTol', 1e-8,'MaxIntervalCount', 20000);
    %  !
    %  !....... Along right side of rectangle
    %  !
    right = quadgk(@myfun,zpoint(npts),zpoint(2*npts-1),...
    'RelTol', 1e-8, 'AbsTol', 1e-8,'MaxIntervalCount', 20000);
    %  !
    %  !....... Along upper side of rectangle
    %  !
    top = quadgk(@myfun,zpoint(2*npts-1),zpoint(3*npts-2),...
    'RelTol', 1e-8, 'AbsTol', 1e-8,'MaxIntervalCount', 20000);
    %  !
    %  !....... Along left side of rectangle
    %  !
    left = quadgk(@myfun,zpoint(3*npts-2),zpoint(4*npts-3),...
    'RelTol', 1e-8, 'AbsTol', 1e-8,'MaxIntervalCount', 20000);
    s(pow) = (lower + right + top + left)/(2*pi*1i);
  end
end
