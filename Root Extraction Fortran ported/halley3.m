function [zroot, bool_converged] = halley3(qpl,qpt, qroot, droot_converged,bool_aitken)
% ***********************************************************************
%
%      This routine implements HALLEY's method, without derivatives,
%      for the computation of zeros of a complex function inside a
%      rectangle in the complex plane.
%
%      The routine implements the iteration formula presented
%      as equation 10-11 of the book
%
%         Iterative Methods for the Solution of Equation
%
%         by Joe Fred Traub (1964)
%
%         Prentice-Hall
%
%      The iteration formula is simply the following
%
%            x(m+1) = x(m) - [ f[x(m)] / theta(m) ]
%
%      where
%
%
%
%
maxiter = 1e2;
%
% -----------------------------------------------------------------------
%
%      I N I T I A L I Z A T I O N
%
% -----------------------------------------------------------------------
%
% ---- Set a ridiculous default value for return
%
zroot = 100;
%
% ---- Set up a small increment either side of the approximate root.
%
%      This is somewhat arbitrary
%
%      The divisor depends on the size of the overall box
%
qdx = (qpt - qpl)/10;
%
% ---- Store the first three estimates of the root
%
z(1)  = qroot;
z(2)  = qroot - qdx;
z(3)  = qroot + qdx;
%
% ---- Store the function value at each estimate of the root
%
zf(1) = FZ(z(1));
zf(2) = FZ(z(2));
zf(3) = FZ(z(3));
%
% ---- Compute the values z(m) - z(m-1) and store in the table at row m
%
delz_m_m1(1)  = 0;
delz_m_m1(2)  = z(2) - z(1);
delz_m_m1(3)  = z(3) - z(2);
%
% ---- Compute the values z(m) - z(m-2) and store in the table at row m
%
delz_m_m2(1)  = 0;
delz_m_m2(2)  = 0;
delz_m_m2(3)  = z(3) - z(1);
%
% ---- Compute the values f( z(m) ) - f( z(m-1) )
%
delf_m_m1(1)  = 0;
delf_m_m1(2)  = zf(2) - zf(1);
delf_m_m1(3)  = zf(3) - zf(2);
%
% ---- Before we proceed, better check if any of our initial
%      estimates are in fact the root
%
for  m = 1 : 3
    if(abs(zf(m)) < droot_converged)
        bool_converged = true;
        zroot          = z(m);
        %
        return;
    end
end
%
% -----------------------------------------------------------------------
%
%      I T E R A T I O N    L O O P
%
% -----------------------------------------------------------------------
%
bool_converged = false;
for iter = 4 : maxiter
    %
    % ....... We are computing point iter --> m+1
    %
    %         So we shall use points m, m-1 and m-2.
    %
    m = iter - 1;
    %
    % ....... So, we now compute D, as defined by Traub, using the
    %         values that exist in our table from previous points.
    %
    %         The first part is f[ z(m), z(m-1) ]
    %
    
    traubsD_first  = delf_m_m1(m) / delz_m_m1(m);
    %
    % ....... The second part is
    %
    %
    %                         f[ z(m), z(m-1), z(m-2) ]
    %           f(z(m-1)) x -----------------------------
    %                             f[ z(m), z(m-1) ]
    %
    %         The fraction here can be expanded out and we use that
    %         form in the calculations
    %
    % ....... The numerator is
    %
    znumer = delz_m_m1(m-1)*delf_m_m1(m) - delz_m_m1(m)*delf_m_m1(m-1);
    %
    % ....... The denominator is
    %
    zdenom = delz_m_m2(m) * delz_m_m1(m-1) * delf_m_m1(m);
    %
    % ....... Now, put the parts together to compute the D defined
    %         by Traub.
    %
    traubsD = traubsD_first - ( zf(m-1)*zdenom/znumer );
    %
    % ....... Now apply equation 10-11, the iteration formula
    %
    z(iter) = z(m) - ( zf(m)/traubsD );
    %
    % ....... Better make sure that this estimate lies in the box
    %
    if(~bool_in_box(qpt,qpl,z(iter)) )
        %    write(6,9900)
        %    write(6,9920) iter,z(iter)
        %     warning('Error: Input outside the box');
    end
    %
    % ....... Ok, now that we have new point let's compute the function
    %         at that point. We'll need this to examine convergence anyway.
    %
    zf(iter) = FZ(z(iter));
    %
    func_abs = abs(zf(iter));
    %
    %    if(bool_debug)then
    %      write(6,2100) iter, z(iter), zf(iter), func_abs
    %    endif
    %
    % ....... Have we converged ?
    %
    if(func_abs < droot_converged)
        %
        %      write(6,2910) droot_converged, iter, z(iter)
        %
        bool_converged = true;
        zroot          = z(iter);
        %
        return
    end
    %
    % ....... Ok, we have not converged yet.
    %
    %         Should we use Aitken acceleration to try again ?
    %
    %         Obviously we do not do this until we have three computed
    %         estimates (maxiter >= 6)
    %
    %         The formula programmed here is from page 266 of the
    %         book by Traub mentioned above, where
    %
    %             x_{i}   = z(iter-2)
    %             x_{i+1} = z(iter-1)
    %             x_{i+2} = z(iter)
    %
    %
    if(bool_aitken)
        zaitken_top       = z(iter-1) - z(iter-2);
        zaitken_top       = zaitken_top*zaitken_top;
        %
        zaitken_bottom    = z(iter) - 2*z(iter-1) + z(iter-2);
        %
        zaitken_estimate  = z(iter-2) - (zaitken_top/zaitken_bottom);
        %
        zf_aitken =  FZ(zaitken_estimate);
        %
        func_abs_ait      = abs(zf_aitken);
        %
        %    if(bool_debug)then
        %    write(6,2200) zaitken_estimate, zf_aitken, func_abs_ait
        % endif
        %
        % ......... Have we converged after Aitken ?
        %
        if(func_abs_ait < droot_converged)
            %
            % write(6,2910) droot_converged, iter, zaitken_estimate
            %
            bool_converged = true;
            zroot          = z(iter);
            %
            return
        end
        %
        % ......... Ok, we have not converged yet but did we get a result
        %           closer to the zero ?
        %
        %           If so, we update the iteration point
        %
        if(func_abs_ait < func_abs)
            z(iter)  = zaitken_estimate;
            zf(iter) = zf_aitken;
        end
    end
    %
    % ....... Ok, for the next iteration we need to fill out
    %         row "iter" of the table.
    %
    delz_m_m1(iter) = z(iter)  - z(iter-1);
    delz_m_m2(iter) = z(iter)  - z(iter-2);
    delf_m_m1(iter) = zf(iter) - zf(iter-1);
    
end
end
