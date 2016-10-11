function [f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last]=DQAGE(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last);
%***BEGIN PROLOGUE  DQAGE
%***PURPOSE  The routine calculates an approximation result to a given
%            definite integral   I = Integral of F over (A,B),
%            hopefully satisfying following claim for accuracy
%            ABS(I-RESLT).LE.MAX(EPSABS,EPSREL*ABS(I)).
%***LIBRARY   SLATEC (QUADPACK)
%***CATEGORY  H2A1A1
%***TYPE      doubleprecision (QAGE-S, DQAGE-D)
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
%        PARAMETERS
%         ON ENTRY
%            F      - doubleprecision
%                     function subprogram defining the integrand
%                     function F(X). The actual name for F needs to be
%                     declared E X T E R N A L in the driver program.
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
%                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
%                     the routine will end with IER = 6.
%
%            KEY    - Integer
%                     Key for choice of local integration rule
%                     A Gauss-Kronrod pair is used with
%                          7 - 15 points if KEY.LT.2,
%                         10 - 21 points if KEY = 2,
%                         15 - 31 points if KEY = 3,
%                         20 - 41 points if KEY = 4,
%                         25 - 51 points if KEY = 5,
%                         30 - 61 points if KEY.GT.5.
%
%            LIMIT  - Integer
%                     Gives an upper bound on the number of subintervals
%                     in the partition of (A,B), LIMIT.GE.1.
%
%         ON RETURN
%            RESULT - doubleprecision
%                     Approximation to the integral
%
%            ABSERR - doubleprecision
%                     Estimate of the modulus of the absolute error,
%                     which should equal or exceed ABS(I-RESULT)
%
%            NEVAL  - Integer
%                     Number of integrand evaluations
%
%            IER    - Integer
%                     IER = 0 Normal and reliable termination of the
%                             routine. It is assumed that the requested
%                             accuracy has been achieved.
%                     IER.GT.0 Abnormal termination of the routine
%                             The estimates for result and error are
%                             less reliable. It is assumed that the
%                             requested accuracy has not been achieved.
%            ERROR MESSAGES
%                     IER = 1 Maximum number of subdivisions allowed
%                             has been achieved. One can allow more
%                             subdivisions by increasing the value
%                             of LIMIT.
%                             However, if this yields no improvement it
%                             is rather advised to analyze the integrand
%                             in order to determine the integration
%                             difficulties. If the position of a local
%                             difficulty can be determined(e.g.
%                             SINGULARITY, DISCONTINUITY within the
%                             interval) one will probably gain from
%                             splitting up the interval at this point
%                             and calling the integrator on the
%                             subranges. If possible, an appropriate
%                             special-purpose integrator should be used
%                             which is designed for handling the type of
%                             difficulty involved.
%                         = 2 The occurrence of roundoff error is
%                             detected, which prevents the requested
%                             tolerance from being achieved.
%                         = 3 Extremely bad integrand behaviour occurs
%                             at some points of the integration
%                             interval.
%                         = 6 The input is invalid, because
%                             (EPSABS.LE.0 and
%                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
%                             RESULT, ABSERR, NEVAL, LAST, RLIST(1) ,
%                             ELIST(1) and IORD(1) are set to zero.
%                             ALIST(1) and BLIST(1) are set to A and B
%                             respectively.
%
%            ALIST   - doubleprecision
%                      Vector of dimension at least LIMIT, the first
%                       LAST  elements of which are the left
%                      end points of the subintervals in the partition
%                      of the given integration range (A,B)
%
%            BLIST   - doubleprecision
%                      Vector of dimension at least LIMIT, the first
%                       LAST  elements of which are the right
%                      end points of the subintervals in the partition
%                      of the given integration range (A,B)
%
%            RLIST   - doubleprecision
%                      Vector of dimension at least LIMIT, the first
%                       LAST  elements of which are the
%                      integral approximations on the subintervals
%
%            ELIST   - doubleprecision
%                      Vector of dimension at least LIMIT, the first
%                       LAST  elements of which are the moduli of the
%                      absolute error estimates on the subintervals
%
%            IORD    - Integer
%                      Vector of dimension at least LIMIT, the first K
%                      elements of which are pointers to the
%                      error estimates over the subintervals,
%                      such that ELIST(IORD(1)), ...,
%                      ELIST(IORD(K)) form a decreasing sequence,
%                      with K = LAST if LAST.LE.(LIMIT/2+2), and
%                      K = LIMIT+1-LAST otherwise
%
%            LAST    - Integer
%                      Number of subintervals actually produced in the
%                      subdivision process
%
%***REFERENCES  (NONE)
%***ROUTINES CALLED  D1MACH, DQK15, DQK21, DQK31, DQK41, DQK51, DQK61,
%                    DQPSRT
%***REVISION HISTORY  (YYMMDD)
%   800101  DATE WRITTEN
%   890531  Changed all specific intrinsics to generic.  (WRB)
%   890831  Modified array declarations.  (WRB)
%   890831  REVISION DATE from Version 3.2
%   891214  Prologue converted to Version 4.0 format.  (BAB)
%***end PROLOGUE  DQAGE
%
persistent a1 a2 area area1 area12 area2 b1 b2 defab1 defab2 defabs epmach errbnd errmax erro12 error1 error2 errsum iroff1 iroff2 k keyf maxerr nrmax resabs uflow ;

if isempty(area), area=0; end;
if isempty(area1), area1=0; end;
if isempty(area12), area12=0; end;
if isempty(area2), area2=0; end;
if isempty(a1), a1=0; end;
if isempty(a2), a2=0; end;
if isempty(b1), b1=0; end;
if isempty(b2), b2=0; end;
if isempty(defabs), defabs=0; end;
if isempty(defab1), defab1=0; end;
if isempty(defab2), defab2=0; end;
if isempty(epmach), epmach=0; end;
if isempty(errbnd), errbnd=0; end;
if isempty(errmax), errmax=0; end;
if isempty(error1), error1=0; end;
if isempty(error2), error2=0; end;
if isempty(erro12), erro12=0; end;
if isempty(errsum), errsum=0; end;
if isempty(resabs), resabs=0; end;
if isempty(uflow), uflow=0; end;
if isempty(iroff1), iroff1=0; end;
if isempty(iroff2), iroff2=0; end;
if isempty(k), k=0; end;
if isempty(keyf), keyf=0; end;
if isempty(maxerr), maxerr=0; end;
if isempty(nrmax), nrmax=0; end;
%
alist_shape=size(alist);alist=reshape(alist,1,[]);
blist_shape=size(blist);blist=reshape(blist,1,[]);
elist_shape=size(elist);elist=reshape(elist,1,[]);
iord_shape=size(iord);iord=reshape(iord,1,[]);
rlist_shape=size(rlist);rlist=reshape(rlist,1,[]);
%
%
%            LIST OF MAJOR VARIABLES
%            -----------------------
%
%           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
%                       CONSIDERED UP TO NOW
%           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
%                       CONSIDERED UP TO NOW
%           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
%                      (ALIST(I),BLIST(I))
%           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
%           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST
%                       ERROR ESTIMATE
%           ERRMAX    - ELIST(MAXERR)
%           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
%           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
%           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
%                       ABS(RESULT))
%           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
%           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
%           LAST      - INDEX FOR SUBDIVISION
%
%
%           MACHINE DEPENDENT CONSTANTS
%           ---------------------------
%
%           EPMACH  IS THE LARGEST RELATIVE SPACING.
%           UFLOW  IS THE SMALLEST POSITIVE MAGNITUDE.
%
%***FIRST EXECUTABLE STATEMENT  DQAGE
[epmach ]=d1mach(4);
[uflow ]=d1mach(1);
%
%           TEST ON VALIDITY OF PARAMETERS
%           ------------------------------
%
ier = 0;
neval = 0;
last = 0;
result = 0.0d+00;
abserr = 0.0d+00;
alist(1) = a;
blist(1) = b;
rlist(1) = 0.0d+00;
elist(1) = 0.0d+00;
iord(1) = 0;
if( epsabs<=0.0d+00 && epsrel<max(0.5d+02.*epmach,0.5d-28) )
ier = 6;
end;
if( ier~=6 )
%
%           FIRST APPROXIMATION TO THE INTEGRAL
%           -----------------------------------
%
keyf = fix(key);
if( key<=0 )
keyf = 1;
end;
if( key>=7 )
keyf = 6;
end;
neval = 0;
if( keyf==1 )
[f,a,b,result,abserr,defabs,resabs]=dqk15(f,a,b,result,abserr,defabs,resabs);
end;
if( keyf==2 )
[f,a,b,result,abserr,defabs,resabs]=dqk21(f,a,b,result,abserr,defabs,resabs);
end;
if( keyf==3 )
[f,a,b,result,abserr,defabs,resabs]=dqk31(f,a,b,result,abserr,defabs,resabs);
end;
if( keyf==4 )
[f,a,b,result,abserr,defabs,resabs]=dqk41(f,a,b,result,abserr,defabs,resabs);
end;
if( keyf==5 )
[f,a,b,result,abserr,defabs,resabs]=dqk51(f,a,b,result,abserr,defabs,resabs);
end;
if( keyf==6 )
[f,a,b,result,abserr,defabs,resabs]=dqk61(f,a,b,result,abserr,defabs,resabs);
end;
last = 1;
rlist(1) = result;
elist(1) = abserr;
iord(1) = 1;
%
%           TEST ON ACCURACY.
%
errbnd = max(epsabs,epsrel.*abs(result));
if( abserr<=0.5d+02.*epmach.*defabs && abserr>errbnd )
ier = 2;
end;
if( limit==1 )
ier = 1;
end;
if( ~(ier~=0 ||(abserr<=errbnd && abserr~=resabs) ||abserr==0.0d+00) )
%
%           INITIALIZATION
%           --------------
%
%
errmax = abserr;
maxerr = 1;
area = result;
errsum = abserr;
nrmax = 1;
iroff1 = 0;
iroff2 = 0;
%
%           MAIN DO-LOOP
%           ------------
%
for last = 2 : limit;
%
%           BISECT THE SUBINTERVAL WITH THE LARGEST ERROR ESTIMATE.
%
a1 = alist(maxerr);
b1 = 0.5d+00.*(alist(maxerr)+blist(maxerr));
a2 = b1;
b2 = blist(maxerr);
if( keyf==1 )
[f,a1,b1,area1,error1,resabs,defab1]=dqk15(f,a1,b1,area1,error1,resabs,defab1);
end;
if( keyf==2 )
[f,a1,b1,area1,error1,resabs,defab1]=dqk21(f,a1,b1,area1,error1,resabs,defab1);
end;
if( keyf==3 )
[f,a1,b1,area1,error1,resabs,defab1]=dqk31(f,a1,b1,area1,error1,resabs,defab1);
end;
if( keyf==4 )
[f,a1,b1,area1,error1,resabs,defab1]=dqk41(f,a1,b1,area1,error1,resabs,defab1);
end;
if( keyf==5 )
[f,a1,b1,area1,error1,resabs,defab1]=dqk51(f,a1,b1,area1,error1,resabs,defab1);
end;
if( keyf==6 )
[f,a1,b1,area1,error1,resabs,defab1]=dqk61(f,a1,b1,area1,error1,resabs,defab1);
end;
if( keyf==1 )
[f,a2,b2,area2,error2,resabs,defab2]=dqk15(f,a2,b2,area2,error2,resabs,defab2);
end;
if( keyf==2 )
[f,a2,b2,area2,error2,resabs,defab2]=dqk21(f,a2,b2,area2,error2,resabs,defab2);
end;
if( keyf==3 )
[f,a2,b2,area2,error2,resabs,defab2]=dqk31(f,a2,b2,area2,error2,resabs,defab2);
end;
if( keyf==4 )
[f,a2,b2,area2,error2,resabs,defab2]=dqk41(f,a2,b2,area2,error2,resabs,defab2);
end;
if( keyf==5 )
[f,a2,b2,area2,error2,resabs,defab2]=dqk51(f,a2,b2,area2,error2,resabs,defab2);
end;
if( keyf==6 )
[f,a2,b2,area2,error2,resabs,defab2]=dqk61(f,a2,b2,area2,error2,resabs,defab2);
end;
%
%           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
%           AND ERROR AND TEST FOR ACCURACY.
%
neval = fix(neval + 1);
area12 = area1 + area2;
erro12 = error1 + error2;
errsum = errsum + erro12 - errmax;
area = area + area12 - rlist(maxerr);
if( defab1~=error1 && defab2~=error2 )
if( abs(rlist(maxerr)-area12)<=0.1d-04.*abs(area12) &&erro12>=0.99d+00.*errmax )
iroff1 = fix(iroff1 + 1);
end;
if( last>10 && erro12>errmax )
iroff2 = fix(iroff2 + 1);
end;
end;
rlist(maxerr) = area1;
rlist(last) = area2;
errbnd = max(epsabs,epsrel.*abs(area));
if( errsum>errbnd )
%
%           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
%
if( iroff1>=6 || iroff2>=20 )
ier = 2;
end;
%
%           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF SUBINTERVALS
%           EQUALS LIMIT.
%
if( last==limit )
ier = 1;
end;
%
%           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
%           AT A POINT OF THE INTEGRATION RANGE.
%
if( max(abs(a1),abs(b2))<=(0.1d+01+0.1d+03.*epmach).*(abs(a2)+0.1d+04.*uflow) )
ier = 3;
end;
end;
%
%           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
%
if( error2>error1 )
alist(maxerr) = a2;
alist(last) = a1;
blist(last) = b1;
rlist(maxerr) = area2;
rlist(last) = area1;
elist(maxerr) = error2;
elist(last) = error1;
else;
alist(last) = a2;
blist(maxerr) = b1;
blist(last) = b2;
elist(maxerr) = error1;
elist(last) = error2;
end;
%
%           CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
%           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
%           WITH THE LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
%
[limit,last,maxerr,errmax,elist,iord,nrmax]=dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax);
% ***JUMP OUT OF DO-LOOP
if( ier~=0 || errsum<=errbnd )
break;
end;
end;
%
%           COMPUTE FINAL RESULT.
%           ---------------------
%
result = 0.0d+00;
for k = 1 : last;
result = result + rlist(k);
end; k = fix(last+1);
abserr = errsum;
end;
if( keyf~=1 )
neval =fix((10.*keyf+1).*(2.*neval+1));
end;
if( keyf==1 )
neval = fix(30.*neval + 15);
end;
end;
alist_shape=zeros(alist_shape);alist_shape(:)=alist(1:numel(alist_shape));alist=alist_shape;
blist_shape=zeros(blist_shape);blist_shape(:)=blist(1:numel(blist_shape));blist=blist_shape;
elist_shape=zeros(elist_shape);elist_shape(:)=elist(1:numel(elist_shape));elist=elist_shape;
iord_shape=zeros(iord_shape);iord_shape(:)=iord(1:numel(iord_shape));iord=iord_shape;
rlist_shape=zeros(rlist_shape);rlist_shape(:)=rlist(1:numel(rlist_shape));rlist=rlist_shape;
end
%DECK DQAG
