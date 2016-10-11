function func_integrand(cside, ctype, fixed_side, ipower, t)
% !***********************************************************************
% !
% !     func_integrand() is passed (as an external reference) to DQAG
% !     and so it called repeatedly from DQAG to evaluate the integrand
% !     at points in the Gauss-Kronod quadrature.
% !
% !     We use DQAG to compute the integral of z**p/f(z) so this routine
% !     computes z**p/f(z), p =0,1,2,3....
% !
% !     It relies on the routine FZ() which computes f(z)
% !
% !     It is useful to remember that we have to integrate around the
% !     box in the complex plane, one side at a time, real part and then
% !     imaginary part. This explains some of the logic seen below.
% !
% !***********************************************************************
% !
% !---- Set default return value
% !
% !
% !---- First of all we have to build the value of z in the complex
% !     plane. This depends on the input "t" and on the side of
% !     the box over which we are integrating.
% !
% !     We also check here to see that cside was initialized properly!
% !
% !     Along RIGHT and LEFT sides  dz = i dy
% !
% !     So we need to multiply integrand by "i". Fixing this here in zmul
% !     avoids another if test later.
% !
% if( (index(cside, 'LOWER') .ne. 0) .or. (index(cside, 'UPPER') .ne. 0) )then
% global cside ctype fixed_side
global F

if (strfind(cside,'LOWER') ||  strfind(cside,'UPPER'))
    
    % !
    z = complex(t, fixed_side);
    % !
    zmul = complex(1, 0);
    % !
elseif (strfind(cside,'RIGHT') ||  strfind(cside,'LEFT'))
    % !
    z = complex(fixed_side, t);
    % !
    zmul = complex(0, 1);
    % !
end

% !
% !---- Make a copy to preserve accross calls into worker routines
% !

ztemp = z;

% !
% !---- We compute the function "f(z)" in the complex plane now
% !

fnz = FZ(z);

% !
% !---- Now compute 1/f(z)
% !

if(ipower ~= 0)
    
    fnz = (ztemp^ipower)/fnz;
    
else
    
    fnz = complex(1, 0)/fnz;
    
end
% !
% !---- Allow for line element having "i" in it.
% !

fnz = fnz*zmul;

% !
% !---- The return value depends side of the rectangle over which
% !     we are integrating.
% !

if(strfind(ctype,'IMAG'))
    
    F = imag(fnz);
    
elseif(strfind(ctype,'REAL'))
    
    F = real(fnz);
    
else
    
    F = 0;
    error('Error. \nInput unspecified');
    
end
end
