function [points, steps, nroots, maxboxes,npts] =  split_contour(point,step,maxroot, npts)
% !***********************************************************************
% !
% !     Given the initial defintion of the rectangular contour, in
% !     point() and Step(), this routine checks to see that we have no
% !     more than "maxroots" roots within the contour. If "maxroots" is
% !     exceeded, then the contour is split internally into smaller
% !     boxes until such times as we have boxes with less than or equal
% !     to "maxroots" roots per box. The rest of the code will then solve
% !     for all of the boxes.
% !
% !***********************************************************************
% % !
% % !---- The function to be integrated
% % !
% % !     Remember that we have to manipulate the real and imaginary
% % !     parts of the integral separately as DQAG is only defined for
% % !     double precision numbers.
% % !
% !
% !---- Values fixed for the duration of the run
% !

% global F
DTWO = 2;
vsmall = 1e-06;
minus1 = -1;
izero = 0;
ione = 1;
s = zeros(1,maxroot*2+1);

% !---- We start by copying the whole contour to the start of the
% !     list of boxes
% !
% !     "nlimit" is a counter that defines the high water marker
% !     for the list of boxes used.
% !
points = point;
steps  = step;
nroots  = minus1;
% !
nlimit = 1;
% !
% !---- We need to set some parameters for the routine DQAG
% !     and then to dynamically allocate the workspace arrays.
% !
% !         "limit" must be .GE. 1 and determines the maximum
% !         number of subintervals in the adaptive partitioning
% !         of the integration interval. It is the dimension of
% !         the array "iwork"
% !
limit = 10000;
% !
% allocate( iwork(1:limit), stat=ierror )
% !
% if(ierror .ne. 0)then
%   write(6,9900)
%   stop 999
% endif
% !
% !.... "lenw" is the dimension of the array "work". It has to be
% !     at least limit*4 in size but we add some head room here just
% !     in case it needs it.
% !
lenw = limit*4 + 100;
% !
% allocate( work(1:lenw), stat=ierror )
% !
% if(ierror .ne. 0)then
%   write(6,9900)
%   stop 999
% endif
% !
% !---- Set the Gauss-Kronod rule required
% !
key = 6;
% !
% !---- Set accuracy requested
% !
epsabs = 1e-07;
epsrel = 1e-07;
% !
% !---- Count the number of roots lying inside the whole contour.
% !
% !     If this is less than or equal to maxroots then we're done
% !     as far as this routine is concerned. Its up to the caller
% !     to handle the case of zero roots !!!
% !
zadaptive = false;
% !
% Call countz
% [nroots(1), s] = countz (points(1,1), Steps(1,1), 'FALSE', npts, maxroot);
[nroots, s] = countz (points, steps, npts, maxroot);
% !
zdone = nroots(1) < maxroot;
% !
% write(6,1500) nroots(1)
% !
% if(zdone)then
%   write(6,1505)
% else
%   write(6,1510)
% endif
% !
% !---- While necessary, split the contour into smaller boxes.
% !
npass = 1;
% !
while(~zdone)
    npass = npass + 1;

    %   write(6,2000) npass, nlimit
    % !
    nbox = nlimit;
    % !
    % !...... Loop over all existing boxes; any that have more
    % !       than "maxroot" get split.
    % !
    for i = 1 : nbox
        if(nroots(i) > maxroot)
            %        write(6,2100) i
            % !
            % !/////////// First IF branch is for split parallel to the IMAGINARY axis
            % !
            if(Steps(1,i) > Steps(2,i))
                t1point(1) = points(1,i);
                t1point(2) = points(2,i);
                t1Step(1)  = Steps(1,i)/DTWO;
                t1Step(2)  = Steps(2,i);
                % !
                t2point(1) = t1point(1) + t1Step(1);
                t2point(2) = t1point(2);
                t2Step(1)  = t1Step(1);
                t2Step(2)  = t1Step(2);
                % !
                % !............. Make a copy of this result, in case we end up moving the contour
                % !              all the way to the right and then have to move it all the way back
                % !              to the left again.
                % !
                t1point_saved(1) = t1point(1);
                t1point_saved(2) = t1point(2);
                t1Step_saved(1)  = t1Step(1);
                t1Step_saved(2)  = t1Step(2);
                % !
                t2point_saved(1) = t2point(1);
                t2point_saved(2) = t2point(2);
                t2Step_saved(1)  = t2Step(1);
                t2Step_saved(2)  = t2Step(2);
                % !
                %          if(zdebug)then
                %            write(6,2102) t1point(1), t1point(2), t1Step(1), t1Step(2), &
                %                          t2point(1), t2point(2), t2Step(1), t2Step(2)
                %          endif
                % !
                % !............. Ok, we'd better check that the new "split contour" that we
                % !              have introduced along the imaginary axis is not too close
                % !              to a root.
                % !
                % !              Contour runs along the line from
                % !
                % !                (t2point(1), t2point(2)) to (t2point(1),      t2point(2)+t2Step(2))
                % !
                % !              We allow movement of ten percent in either direction.
                % !
                zfinished_adjust = false;
                nadjusts         = 0;
                zdelta_Step      = complex(t1Step(1)/(10*real(100)), 0);
                % !
                %  while(.not.zfinished_adjust)
                while(~zfinished_adjust)
                    zstart  = complex(t2point(1), t2point(2));
                    zfinish = complex(t2point(1), t2point(2) + t2Step(2));
                    % !
                    % !................ Debug printout of the contour along which we shall be
                    % !                 performing the integration.
                    % !
                    %             if(zdebug)then
                    %               write(6,2120) nadjusts, zstart, zfinish
                    %             endif
                    % !
                    % !................ Fot this iteration, we set a flag to denote whether we
                    % !                 need to adjust the contour of not.
                    % !
                    zclose_to_root = false;
                    % !
                    % !................ Now check if we need to adjust the contour.
                    % !
                    % !                 We use DQAG to check convergence of the integral  1/f(z)
                    % !
                    % !                 Start with the REAL part
                    % !
                    np = 0;
                    % !
                    % ier = ifunc_integrand_setup('RIGHT', t1point(1)+t1Step(1), 'REAL', np);
                    ier = ifunc_integrand_setup('RIGHT',t1point(1)+t1Step(1),'REAL', np);
                    % !
                    % call DQAG(func_integrand, t1point(2), t1point(2)+t1Step(2),      &
                    %           epsabs, epsrel, key,                                   &
                    %           result_right_real, abserr,  neval,  ier,               &
                    %           limit, lenw, last, iwork, work)
                    [result_right_real, abserr,  neval,  ier,       ...
                        limit, lenw, last, iwork, work] =  DQAG(func_integrand, t1point(2), t1point(2) + t1Step(2),  ...
                        epsabs, epsrel, key);
                    %  [f,a,b,epsabs,epsrel,key,result_right_real,abserr,neval,ier,limit,lenw,last,iwork,work]=dqag(func_integrand,t1point(2),t1point(2),epsabs,epsrel,key,result,abserr,neval,ier,limit,lenw,last,iwork,work);
                    % !
                    %             if(zdebug)then
                    %               write(6,3010)  ier, 'REAL', abserr, neval
                    %               if(ier .eq. 1) write(6,3001)
                    %               if(ier .eq. 2) write(6,3002)
                    %               if(ier .eq. 3) write(6,3003)
                    %               if(ier .eq. 6) write(6,3006)
                    %             endif
                    % !
                    zclose_to_root = ier ~= 0;
                    % !
                    % !................ Now the IMAG part
                    % !
                    np = 0;
                    % !
                    ier = ifunc_integrand_setup('RIGHT', t1point(1)+t1Step(1), 'IMAG', np);
                    % !
                    % call DQAG(func_integrand, t1point(2), t1point(2)+t1Step(2),      &
                    %           epsabs, epsrel, key,                                   &
                    %           result_right_imag, abserr,  neval,  ier,               &
                    %           limit, lenw, last, iwork, work)

                    [result_right_imag, abserr,  neval,  ier,       ...
                        limit, lenw, last, iwork, work] =  DQAG(func_integrand, t1point(2), t1point(2)+t1Step(2),  ...
                        epsabs, epsrel, key);
                    % !
                    %             if(zdebug)then
                    %               write(6,3010)  ier, 'IMAG', abserr, neval
                    %               if(ier .eq. 1) write(6,3001)
                    %               if(ier .eq. 2) write(6,3002)
                    %               if(ier .eq. 3) write(6,3003)
                    %               if(ier .eq. 6) write(6,3006)
                    %             endif
                    % !
                    zclose_to_root = ier ~= 0;
                    % !
                    % !................ Finally, scan along the line for very small values of ABS(f(z))
                    % !
                    % !                 If there is a root on the line, DQAG does not report an error!
                    % !
                    % !                 Need to scan very tightly here indeed.
                    % !
                    zdelta = (zfinish - zstart)/complex(real(3*npts-1),0);
                    % !
                    jpoint = 1;
                    % !
                    while((jpoint <= 3*npts-1) && ~zclose_to_root)
                        zpoint = zstart + complex(real(jpoint),0)*zdelta;
                        % !
                        zfunc = FZ(zpoint);
                        % !
                        if(abs(zfunc) < 1e-05)
                            %  if(zdebug)then
                            %     write(6,2121) jpoint, zpoint, ABS(zfunc)
                            %  endif
                            zclose_to_root = true;
                        end
                        % !
                        jpoint = jpoint + 1;
                    end
                    % !
                    % !................ Did either integral trip an error condition perhaps ???
                    % !
                    %             if(zdebug)then
                    %               write(6,2125) zclose_to_root
                    %             endif
                    % !
                    % !................ Ok, if we do not need to adjust the contour then set the
                    % !                 flag to reflect this, otherwise adjust and try again.
                    % !
                    % !                  We go to the right and then to the left
                    % !
                    if(~zclose_to_root)
                        zfinished_adjust = true;
                    else
                        if(nadjusts < 100)
                            nadjusts   = nadjusts + 1;
                            t1Step(1)  = t1Step(1)  + zdelta_Step;
                            t2point(1) = t1point(1) + t1Step(1);
                            t2Step(1)  = t2Step(1)  - zdelta_Step;
                        elseif(nadjusts == 100)
                            nadjusts   = nadjusts + 1;
                            % !
                            t1point(1) = t1point_saved(1);
                            t1point(2) = t1point_saved(2);
                            t1Step(1)  = t1Step_saved(1);
                            t1Step(2)  = t1Step_saved(2);
                            % !
                            t2point(1) = t2point_saved(1);
                            t2point(2) = t2point_saved(2);
                            t2Step(1)  = t2Step_saved(1);
                            t2Step(2)  = t2Step_saved(2);
                            % !
                            t1Step(1)  = t1Step(1)  - zdelta_Step;
                            t2point(1) = t1point(1) + t1Step(1);
                            t2Step(1)  = t2Step(1)  + zdelta_Step;
                        elseif(nadjusts > 100 && nadjusts <= 200)
                            nadjusts   = nadjusts + 1;
                            t1Step(1)  = t1Step(1)  - zdelta_Step;
                            t2point(1) = t1point(1) + t1Step(1);
                            t2Step(1)  = t2Step(1)  + zdelta_Step;
                        elseif(nadjusts > 200)
                            % write(6,9900)
                            % write(6,9910)
                            % write(6,9915) npass, nbox, nadjusts
                            error('Error. \nInput unspecified');
                        else
                            % write(6,9900)
                            % write(6,9999)
                            % write(6,9915) npass, nbox, nadjusts
                            error('Error. \nInput unspecified');
                        end
                    end
                end
                % !
                % !/////////// Following code is for split parallel to the REAL axis.
                % !
            else
                t1point(1) = points(1,i);
                t1point(2) = points(2,i);
                t1Step(1)  = Steps(1,i);
                t1Step(2)  = Steps(2,i)/DTWO;
                t2point(1) = t1point(1);
                t2point(2) = t1point(2) + t1Step(2);
                t2Step(1)  = t1Step(1);
                t2Step(2)  = t1Step(2);
                % !
                % !............. Make a copy of this result, in case we end up moving the contour
                % !              all the way up and then have to move it all the way back down again.
                % !
                t1point_saved(1) = t1point(1);
                t1point_saved(2) = t1point(2);
                t1Step_saved(1)  = t1Step(1);
                t1Step_saved(2)  = t1Step(2);
                % !
                t2point_saved(1) = t2point(1);
                t2point_saved(2) = t2point(2);
                t2Step_saved(1)  = t2Step(1);
                t2Step_saved(2)  = t2Step(2);
                % !
                %  if(zdebug)then
                %    write(6,2103) t1point(1), t1point(2), t1Step(1), t1Step(2), &
                %                  t2point(1), t2point(2), t2Step(1), t2Step(2)
                %  endif
                % !
                % !............. Ok, we'd better check that the new "split contour" that we
                % !              have introduced along the imaginary axis is not too close
                % !              to a root.
                % !
                % !              Contour runs along the line from
                % !
                % !                (t2point(1), t2point(2)) to (t2point(1) + t2Step(1), t2point(2))
                % !
                % !              We allow adjustment of ten percent in either direction in 100 moves.
                % !
                zfinished_adjust = false;
                nadjusts         = 0;
                zdelta_Step      = complex(0, t2Step(2)/(10*real(100)));
                % !
                while(~zfinished_adjust)
                    zstart  = complex(t2point(1), t2point(2));
                    zfinish = complex(t2point(1) + t2Step(1), t2point(2));
                    % !
                    % !................ Debug printout of the contour along which we shall be
                    % !                 performing the integration.
                    % !
                    %             if(zdebug)then
                    %               write(6,2120) nadjusts, zstart, zfinish
                    %             endif
                    % !
                    % !................ For this iteration, we set a flag to denote whether we
                    % !                 need to adjust the contour of not.
                    % !
                    zclose_to_root =false;
                    % !
                    % !................ Now check if we need to adjust the contour.
                    % !
                    % !                 We use DQAG to check convergence of the integral  1/f(z)
                    % !
                    % !                 Start with the REAL part
                    % !
                    np = 0;
                    % !
                    ier = ifunc_integrand_setup('LOWER', t2point(2), 'REAL', np);
                    % !
                    % call DQAG(func_integrand, t2point(1), t2point(1)+t2Step(1),      &
                    %           epsabs, epsrel, key,                                   &
                    %           result_lower_real, abserr,  neval,  ier,               &
                    %           limit, lenw, last, iwork, work)

                    [result_lower_real, abserr,  neval,  ier,       ...
                        limit, lenw, last, iwork, work] =  DQAG(func_integrand, t2point(1), t2point(1)+t2Step(1),  ...
                        epsabs, epsrel, key);
                    % !
                    %             if(zdebug)then
                    %               write(6,3010)  ier, 'REAL', abserr, neval
                    %               if(ier .eq. 1) write(6,3001)
                    %               if(ier .eq. 2) write(6,3002)
                    %               if(ier .eq. 3) write(6,3003)
                    %               if(ier .eq. 6) write(6,3006)
                    %             endif
                    % !
                    zclose_to_root = ier ~= 0;
                    % !
                    % !................ IMAG part
                    % !
                    np = 0;
                    % !
                    ier = ifunc_integrand_setup('LOWER', t2point(2), 'IMAG', np);
                    % !
                    %             call DQAG(func_integrand, t2point(1), t2point(1)+t2Step(1),      &
                    %                       epsabs, epsrel, key,                                   &
                    %                       result_lower_imag, abserr,  neval,  ier,               &
                    %                       limit, lenw, last, iwork, work)
                    % !
                    [result_lower_imag, abserr,  neval,  ier,       ...
                        limit, lenw, last, iwork, work] =  DQAG(func_integrand, t2point(1), t2point(1)+t2Step(1),  ...
                        epsabs, epsrel, key);
                    %             if(zdebug)then
                    %               write(6,3010)  ier, 'IMAG', abserr, neval
                    %               if(ier .eq. 1) write(6,3001)
                    %               if(ier .eq. 2) write(6,3002)
                    %               if(ier .eq. 3) write(6,3003)
                    %               if(ier .eq. 6) write(6,3006)
                    %             endif
                    % !
                    zclose_to_root = ier ~= 0;
                    % !
                    % !................ Finally, scan along the line for very small values of ABS(f(z))
                    % !
                    % !                 If there is a root on the line, DQAG does not report an error!
                    % !
                    % !                 Need to scan very tightly here indeed.
                    % !
                    zdelta = (zfinish - zstart)/complex(real(3*npts-1),0);
                    % !/
                    jpoint = 1;
                    % !
                    while((jpoint <= 3*npts-1) && ~zclose_to_root)
                        zpoint = zstart + complex(real(jpoint),0)*zdelta;
                        % !
                        zfunc = FZ(zpoint);
                        % !
                        if(abs(zfunc) < 1e-05)
                            %  if(zdebug)then
                            %     write(6,2121) jpoint, zpoint, ABS(zfunc)
                            %  endif
                            zclose_to_root = true;
                        end
                        % !
                        jpoint = jpoint + 1;
                    end
                    % !
                    % !................ Did either integral trip an error condition perhaps ???
                    % !
                    % if(zdebug)then
                    %   write(6,2125) zclose_to_root
                    % endif
                    % !
                    % !................ Ok, if we do not need to adjust the contour then set the
                    % !                 flag to reflect this, otherwise adjust and try again.
                    % !
                    % !                  We go upwards and then downwards
                    % !
                    if(~zclose_to_root)
                        zfinished_adjust = true;
                    else
                        if(nadjusts < 100)
                            nadjusts   = nadjusts + 1;
                            t1Step(2)  = t1Step(2)  + imag(zdelta_Step);
                            t2point(2) = t1point(2) + t1Step(2);
                            t2Step(2)  = t2Step(2)  - imag(zdelta_Step);
                        elseif(nadjusts == 100)
                            nadjusts   = nadjusts + 1;
                            % !
                            t1point(1) = t1point_saved(1);
                            t1point(2) = t1point_saved(2);
                            t1Step(1)  = t1Step_saved(1);
                            t1Step(2)  = t1Step_saved(2);
                            % !
                            t2point(1) = t2point_saved(1);
                            t2point(2) = t2point_saved(2);
                            t2Step(1)  = t2Step_saved(1);
                            t2Step(2)  = t2Step_saved(2);
                            % !
                            t1Step(2)  = t1Step(2)  - imag(zdelta_Step);
                            t2point(2) = t1point(2) + t1Step(2);
                            t2Step(2)  = t2Step(2)  + imag(zdelta_Step);
                        elseif(nadjusts < 100 && nadjusts <= 200)
                            nadjusts   = nadjusts + 1;
                            t1Step(2)  = t1Step(2)  - imag(zdelta_Step);
                            t2point(2) = t1point(2) + t1Step(2);
                            t2Step(2)  = t2Step(2)  + imag(zdelta_Step);
                        elseif(nadjusts > 200)
                            % write(6,9900)
                            % write(6,9910)
                            % write(6,9915) npass, nbox, nadjusts
                            error('Error. \nInput unspecified');
                        else
                            % write(6,9900)
                            % write(6,9999)
                            % write(6,9915) npass, nbox, nadjusts
                            error('Error. \nInput unspecified');
                        end
                    end
                end
                % !
                % !////////// End of SPLIT of boxes
                % !
            end
            % !
            % !........ Re compute roots in the two boxes
            % !
            zadaptive = false;
            % !
            %     call countz(t1point, t1Step, zadaptive, npts, maxroot, nroots1, s)
            %     call countz(t2point, t2Step, zadaptive, npts, maxroot, nroots2, s)
            [nroots1, s] = countz (t1point, t1Step, zadapative, npts, maxroots);
            [nroots2, s] = countz (t2point, t2Step, zadapative, npts, maxroots);
            % !
            %     if(zdebug)then
            %       write(6,2200) 'After contour readjustment, Sub-Box 1 of box ',i
            %       write(6,2210) t1point, t1Step, nroots1
            %       write(6,2200) 'After contour readjustment, Sub-Box 2 of box ',i
            %       write(6,2210) t2point, t2Step, nroots2
            %     endif
            % !
            % !........ Add new box to end of list, remove existing box
            % !
            if((nroots1 ~= 0) && (nroots2 == 0))
                points(1,i) = t1point(1);
                points(2,i) = t1point(2);
                Steps(1,i)  = t1Step(1);
                Steps(2,i)  = t1Step(2);
                nroots(i)   = nroots1;
            elseif((nroots1 == 0) && (nroots2 ~= 0))
                points(1,i) = t2point(1);
                points(2,i) = t2point(2);
                Steps(1,i)  = t2Step(1);
                Steps(2,i)  = t2Step(2);
                nroots(i)   = nroots2;
            elseif((nroots1 ~= 0) && (nroots2 ~= 0))
                if(nlimit+1 <= maxboxes)
                    nlimit = nlimit + 1;
                    % !
                    points(1,i) = t1point(1);
                    points(2,i) = t1point(2);
                    Steps(1,i)  = t1Step(1);
                    Steps(2,i)  = t1Step(2);
                    nroots(i)   = nroots1;
                    % !
                    points(1,nlimit) = t2point(1);
                    points(2,nlimit) = t2point(2);
                    Steps(1,nlimit)  = t2Step(1);
                    Steps(2,nlimit)  = t2Step(2);
                    nroots(nlimit)   = nroots2;
                else
                    error('Error. \nInput unspecified');
                end
            else
                % !
                % !--- Both boxes have zero roots - should not happen
                % !
                % write(6,*) ' ***** Error in: split_contour()'
                % write(6,*) ' both new boxes have zero roots - should not happen'
                error('Error. ***** Error in: split_contour()\nboth new boxes have zero roots - should not happen');
            end
        end
        % !
        % 100   continue
        % !
        % !------ Check if all boxes satisfy root limit.
        % !
        zdone = true;
        % !
        for  i = 1 : nlimit
            if(nroots(i) > maxroot)
                zdone = false;
            end
        end

    end

end
