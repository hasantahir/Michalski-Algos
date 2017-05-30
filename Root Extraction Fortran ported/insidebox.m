function check = insidebox(qpt,qpl,z)
% ***********************************************************************
%
%      This routine determines whether the value z lies inside the
%      rectangular box defined by points qpt, qpl in the complex
%      plane
%
%                                                              qpt
%          +----------------------------------------------------+
%          |                                                    |
%          |                                                    |
%          |                                                    |
%          |                                                    |
%          |                         +                          |
%          |                         z                          |
%          |                                                    |
%          |                                                    |
%          |                                                    |
%          +----------------------------------------------------+
%         qpl
%
%
% ***********************************************************************
%
% ---- Perform test now
%
check = false;

if ( real(z) > real(qpl) && real(z) < real(qpt) && ...
        imag(z) > imag(qpl) && imag(z) < imag(qpt))
    %
    % ---- We get here if we are inside the rectangle
    %
    check = true;
end


end
