function [ifunc_integrand_setup, fixed_side, ipower, cside, ctype] = ifunc_integrand_setup(cside_in,fixed_side_in,ctype_in, ipower_in)
% !=======================================================================
% !
% !     This entry point is called before the routine DQAG is invoked
% !     in order to store necessary details of the side of the box along
% !     which the integrand is being performed. An alternative would
% !     have been a common block but that is more prone to error.
% !
% !=======================================================================
% !

% !
% !---- Let's assume that the input data is correct and only
% !     reset for an error flag then. Note that we invalidate
% !     "cside" in local storage too - it will be used as an
% !     error check.
% !
ifunc_integrand_setup = 1;
% !
% !---- Copy the fixed side to local storage
% !
fixed_side = fixed_side_in;
% !
% !---- Copy the power to local storage
% !
ipower = ipower_in;
% !
% !---- Ok, so check that the "side" is properly defined
% !
if (~isempty(strfind(cside_in,'LOWER')) ~= 0)

    cside = 'LOWER';

elseif (~isempty(strfind(cside_in,'UPPER')) ~= 0)

    cside = 'UPPER';

elseif (~isempty(strfind(cside_in,'RIGHT')) ~= 0)

    cside = 'RIGHT';

elseif(~isempty(strfind(cside_in,'LEFT')) ~= 0)

    cside = 'LEFT';

else

    cside = 'ERROR';
    ifunc_integrand_setup = 0;
    return;

end

% !
% !---- Ok, so check that the "type" is properly defined
% !
if(~isempty(strfind(ctype_in,'IMAG')) ~= 0)

    ctype = 'IMAG';

elseif(~isempty(strfind(ctype_in,'REAL')) ~= 0)

    ctype = 'REAL';

else

    ctype = 'ERROR';
    ifunc_integrand_setup = 0;
    return;

end

end
