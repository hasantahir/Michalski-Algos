function [f,a,b,epsabs,epsrel,key,result,abserr,neval,ier,limit,lenw,last,iwork,work]=dqag(f,a,b,epsabs,epsrel,key,result,abserr,neval,ier,limit,lenw,last,iwork,work);
%***BEGIN PROLOGUE  DQAG
%***PURPOSE  The routine calculates an approximation result to a given
%            definite integral I = integral of F over (A,B),
%            hopefully satisfying following claim for accuracy
%            ABS(I-RESULT)LE.MAX(EPSABS,EPSREL*ABS(I)).
%***LIBRARY   SLATEC (QUADPACK)
%***CATEGORY  H2A1A1
%***TYPE      doubleprecision (QAG-S, DQAG-D)
%***KEYWORDS  AUTOMATIC INTEGRATOR, GAUSS-KRONROD RULES,
%             GENERAL-PURPOSE, GLOBALLY ADAPTIVE, INTEGRAND EXAMINATOR,
%             QUADPACK, QUADRATURE
%***AUTHOR  Piessens, Robert
%             Applied Mathematics and Programming Division
%             K. U. Leuven
%           de Doncker, Elise
%             Applied Mathematics and Programming Division
%             K. U. Leuven
%***DESCRIPTION
%
%        Computation of a definite integral
%        Standard fortran subroutine
%        doubleprecision version
%
%            F      - doubleprecision
%                     function subprogram defining the integrand
%                     function F(X). The actual name for F needs to be
%                     Declared E X T E R N A L in the driver program.
%
%            A      - doubleprecision
%                     Lower limit of integration
%
%            B      - doubleprecision
%                     Upper limit of integration
%
%            EPSABS - doubleprecision
%                     Absolute accuracy requested
%            EPSREL - doubleprecision
%                     Relative accuracy requested
%                     If  EPSABS.LE.0
%                     And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
%                     The routine will end with IER = 6.
%
%            KEY    - Integer
%                     Key for choice of local integration rule
%                     A GAUSS-KRONROD PAIR is used with
%                       7 - 15 POINTS If KEY.LT.2,
%                      10 - 21 POINTS If KEY = 2,
%                      15 - 31 POINTS If KEY = 3,
%                      20 - 41 POINTS If KEY = 4,
%                      25 - 51 POINTS If KEY = 5,
%                      30 - 61 POINTS If KEY.GT.5.
%
%         ON RETURN
%            RESULT - doubleprecision
%                     Approximation to the integral
%
%            ABSERR - doubleprecision
%                     Estimate of the modulus of the absolute error,
%                     Which should EQUAL or EXCEED ABS(I-RESULT)
%
%            NEVAL  - Integer
%                     Number of integrand evaluations
%
%            IER    - Integer
%                     IER = 0 Normal and reliable termination of the
%                             routine. It is assumed that the requested
%                             accuracy has been achieved.
%                     IER.GT.0 Abnormal termination of the routine
%                             The estimates for RESULT and ERROR are
%                             Less reliable. It is assumed that the
%                             requested accuracy has not been achieved.
%                      ERROR MESSAGES
%                     IER = 1 Maximum number of subdivisions allowed
%                             has been achieved. One can allow more
%                             subdivisions by increasing the value of
%                             LIMIT (and taking the according dimension
%                             adjustments into account). HOWEVER, If
%                             this yield no improvement it is advised
%                             to analyze the integrand in order to
%                             determine the integration difficulties.
%                             If the position of a local difficulty can
%                             be determined (I.E. SINGULARITY,
%                             DISCONTINUITY WITHIN THE INTERVAL) One
%                             will probably gain from splitting up the
%                             interval at this point and calling the
%                             INTEGRATOR on the SUBRANGES. If possible,
%                             AN APPROPRIATE SPECIAL-PURPOSE INTEGRATOR
%                             should be used which is designed for
%                             handling the type of difficulty involved.
%                         = 2 The occurrence of roundoff error is
%                             detected, which prevents the requested
%                             tolerance from being achieved.
%                         = 3 Extremely bad integrand behaviour occurs
%                             at some points of the integration
%                             interval.
%                         = 6 The input is invalid, because
%                             (EPSABS.LE.0 AND
%                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
%                             OR LIMIT.LT.1 OR LENW.LT.LIMIT*4.
%                             RESULT, ABSERR, NEVAL, LAST are set
%                             to zero.
%                             EXCEPT when LENW is invalid, IWORK(1),
%                             WORK(LIMIT*2+1) and WORK(LIMIT*3+1) are
%                             set to zero, WORK(1) is set to A and
%                             WORK(LIMIT+1) to B.
%
%         DIMENSIONING PARAMETERS
%            LIMIT - Integer
%                    Dimensioning parameter for IWORK
%                    Limit determines the maximum number of subintervals
%                    in the partition of the given integration interval
%                    (A,B), LIMIT.GE.1.
%                    If LIMIT.LT.1, the routine will end with IER = 6.
%
%            LENW  - Integer
%                    Dimensioning parameter for work
%                    LENW must be at least LIMIT*4.
%                    IF LENW.LT.LIMIT*4, the routine will end with
%                    IER = 6.
%
%            LAST  - Integer
%                    On return, LAST equals the number of subintervals
%                    produced in the subdivision process, which
%                    determines the number of significant elements
%                    actually in the WORK ARRAYS.
%
%         WORK ARRAYS
%            IWORK - Integer
%                    Vector of dimension at least limit, the first K
%                    elements of which contain pointers to the error
%                    estimates over the subintervals, such that
%                    WORK(LIMIT*3+IWORK(1)),... , WORK(LIMIT*3+IWORK(K))
%                    form a decreasing sequence with K = LAST If
%                    LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST otherwise
%
%            WORK  - doubleprecision
%                    Vector of dimension at least LENW
%                    on return
%                    WORK(1), ..., WORK(LAST) contain the left end
%                    points of the subintervals in the partition of
%                     (A,B),
%                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain the
%                     right end points,
%                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
%                     the integral approximations over the subintervals,
%                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST) contain
%                     the error estimates.
%
%***REFERENCES  (NONE)
%***ROUTINES CALLED  DQAGE, XERMSG
%***REVISION HISTORY  (YYMMDD)
%   800101  DATE WRITTEN
%   890831  Modified array declarations.  (WRB)
%   890831  REVISION DATE from Version 3.2
%   891214  Prologue converted to Version 4.0 format.  (BAB)
%   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
%***end PROLOGUE  DQAG
persistent l1 l2 l3 lvl ;

if isempty(lvl), lvl=0; end;
if isempty(l1), l1=0; end;
if isempty(l2), l2=0; end;
if isempty(l3), l3=0; end;
%
iwork_shape=size(iwork);iwork=reshape(iwork,1,[]);
work_shape=size(work);work=reshape(work,1,[]);
%
%***FIRST EXECUTABLE STATEMENT  DQAG
ier = 6;
neval = 0;
last = 0;
result = 0.0d+00;
abserr = 0.0d+00;
if( limit>=1 && lenw>=limit.*4 )
%
%        PREPARE CALL FOR DQAGE.
%
l1 = fix(limit + 1);
l2 = fix(limit + l1);
l3 = fix(limit + l2);
%
[f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,ier,dumvar12,dumvar13,dumvar14,dumvar15,iwork,last]=dqage(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,ier,work(sub2ind(size(work),max(1,1)):end),work(sub2ind(size(work),max(l1,1)):end),work(sub2ind(size(work),max(l2,1)):end),work(sub2ind(size(work),max(l3,1)):end),iwork,last);   dumvar12i=find((work(sub2ind(size(work),max(1,1)):end))~=(dumvar12));dumvar13i=find((work(sub2ind(size(work),max(l1,1)):end))~=(dumvar13));dumvar14i=find((work(sub2ind(size(work),max(l2,1)):end))~=(dumvar14));dumvar15i=find((work(sub2ind(size(work),max(l3,1)):end))~=(dumvar15));   work(1-1+dumvar12i)=dumvar12(dumvar12i); work(l1-1+dumvar13i)=dumvar13(dumvar13i); work(l2-1+dumvar14i)=dumvar14(dumvar14i); work(l3-1+dumvar15i)=dumvar15(dumvar15i);
%
%        CALL ERROR HANDLER IF NECESSARY.
%
lvl = 0;
end;
%
if( ier==6 )
lvl = 1;
end;
if( ier~=0 )
[dumvar1,dumvar2,dumvar3,ier,lvl]=xermsg('SLATEC','DQAG','ABNORMAL RETURN',ier,lvl);
end;
iwork_shape=zeros(iwork_shape);iwork_shape(:)=iwork(1:numel(iwork_shape));iwork=iwork_shape;
work_shape=zeros(work_shape);work_shape(:)=work(1:numel(work_shape));work=work_shape;
end
%DECK DQAGIE
