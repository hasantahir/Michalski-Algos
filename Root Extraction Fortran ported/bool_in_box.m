function bool_rc = bool_in_box(qpt,qpl,z)
% !***********************************************************************
% !
% !     This routine determines whether the value z lies inside the
% !     rectangular box defined by points qpt, qpl in the complex
% !     plane
% !
% !                                                             qpt
% !         +----------------------------------------------------+
% !         |                                                    |
% !         |                                                    |
% !         |                                                    |
% !         |                                                    |
% !         |                         +                          |
% !         |                         z                          |
% !         |                                                    |
% !         |                                                    |
% !         |                                                    |
% !         +----------------------------------------------------+
% !        qpl
% !
% !
% !***********************************************************************
% !
% !---- Default return value
% !
% !     Let's assume the worst - we are outside the box.
% !
bool_rc = false;
% !
% !---- Perform test now
% !
if( real(z) > real(qpt) )
    return;
end
% !
if( real(z) < real(qpl) )
    return;
end
% !
if( imag(z) > imag(qpt) )
    return;
end
% !
if( imag(z) < imag(qpl) )
    return;
end
% !
% !---- We get here if we are inside the rectangle
% !
bool_rc = true;
% !
% !---- Return point
% !
end
