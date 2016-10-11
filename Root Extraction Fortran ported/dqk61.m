function [f,a,b,result,abserr,resabs,resasc]=dqk61(f,a,b,result,abserr,resabs,resasc);
%***BEGIN PROLOGUE  DQK61
%***PURPOSE  To compute I = Integral of F over (A,B) with error
%                           estimate
%                       J = Integral of ABS(F) over (A,B)
%***LIBRARY   SLATEC (QUADPACK)
%***CATEGORY  H2A1A2
%***TYPE      doubleprecision (QK61-S, DQK61-D)
%***KEYWORDS  61-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
%***AUTHOR  Piessens, Robert
%             Applied Mathematics and Programming Division
%             K. U. Leuven
%           de Doncker, Elise
%             Applied Mathematics and Programming Division
%             K. U. Leuven
%***DESCRIPTION
%
%        Integration rule
%        Standard fortran subroutine
%        doubleprecision version
%
%
%        PARAMETERS
%         ON ENTRY
%           F      - doubleprecision
%                    function subprogram defining the integrand
%                    function F(X). The actual name for F needs to be
%                    declared E X T E R N A L in the calling program.
%
%           A      - doubleprecision
%                    Lower limit of integration
%
%           B      - doubleprecision
%                    Upper limit of integration
%
%         ON RETURN
%           RESULT - doubleprecision
%                    Approximation to the integral I
%                    RESULT is computed by applying the 61-point
%                    Kronrod rule (RESK) obtained by optimal addition of
%                    abscissae to the 30-point Gauss rule (RESG).
%
%           ABSERR - doubleprecision
%                    Estimate of the modulus of the absolute error,
%                    which should equal or exceed ABS(I-RESULT)
%
%           RESABS - doubleprecision
%                    Approximation to the integral J
%
%           RESASC - doubleprecision
%                    Approximation to the integral of ABS(F-I/(B-A))
%
%***REFERENCES  (NONE)
%***ROUTINES CALLED  D1MACH
%***REVISION HISTORY  (YYMMDD)
%   800101  DATE WRITTEN
%   890531  Changed all specific intrinsics to generic.  (WRB)
%   890531  REVISION DATE from Version 3.2
%   891214  Prologue converted to Version 4.0 format.  (BAB)
%***end PROLOGUE  DQK61
%
persistent centr dabsc dhlgth epmach fc firstCall fsum fv1 fv2 fval1 fval2 hlgth j jtw jtwm1 resg resk reskh uflow wg wgk xgk ; if isempty(firstCall),firstCall=1;end;

if isempty(dabsc), dabsc=0; end;
if isempty(centr), centr=0; end;
if isempty(dhlgth), dhlgth=0; end;
if isempty(epmach), epmach=0; end;
if isempty(fc), fc=0; end;
if isempty(fsum), fsum=0; end;
if isempty(fval1), fval1=0; end;
if isempty(fval2), fval2=0; end;
if isempty(fv1), fv1=zeros(1,30); end;
if isempty(fv2), fv2=zeros(1,30); end;
if isempty(hlgth), hlgth=0; end;
if isempty(resg), resg=0; end;
if isempty(resk), resk=0; end;
if isempty(reskh), reskh=0; end;
if isempty(uflow), uflow=0; end;
if isempty(wg), wg=zeros(1,15); end;
if isempty(wgk), wgk=zeros(1,31); end;
if isempty(xgk), xgk=zeros(1,31); end;
if isempty(j), j=0; end;
if isempty(jtw), jtw=0; end;
if isempty(jtwm1), jtwm1=0; end;
%
%
%           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE
%           INTERVAL (-1,1). BECAUSE OF SYMMETRY ONLY THE POSITIVE
%           ABSCISSAE AND THEIR CORRESPONDING WEIGHTS ARE GIVEN.
%
%           XGK   - ABSCISSAE OF THE 61-POINT KRONROD RULE
%                   XGK(2), XGK(4)  ... ABSCISSAE OF THE 30-POINT
%                   GAUSS RULE
%                   XGK(1), XGK(3)  ... OPTIMALLY ADDED ABSCISSAE
%                   TO THE 30-POINT GAUSS RULE
%
%           WGK   - WEIGHTS OF THE 61-POINT KRONROD RULE
%
%           WG    - WEIGHTS OF THE 30-POINT GAUSS RULE
%
%
% GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
% AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
% BELL LABS, NOV. 1981.
%
if firstCall,   wg(1)=[0.007968192496166605615465883474674d0];  end;
if firstCall,   wg(2)=[0.018466468311090959142302131912047d0];  end;
if firstCall,   wg(3)=[0.028784707883323369349719179611292d0];  end;
if firstCall,   wg(4)=[0.038799192569627049596801936446348d0];  end;
if firstCall,   wg(5)=[0.048402672830594052902938140422808d0];  end;
if firstCall,   wg(6)=[0.057493156217619066481721689402056d0];  end;
if firstCall,   wg(7)=[0.065974229882180495128128515115962d0];  end;
if firstCall,   wg(8)=[0.073755974737705206268243850022191d0];  end;
if firstCall,   wg(9)=[0.080755895229420215354694938460530d0];  end;
if firstCall,   wg(10)=[0.086899787201082979802387530715126d0];  end;
if firstCall,   wg(11)=[0.092122522237786128717632707087619d0];  end;
if firstCall,   wg(12)=[0.096368737174644259639468626351810d0];  end;
if firstCall,   wg(13)=[0.099593420586795267062780282103569d0];  end;
if firstCall,   wg(14)=[0.101762389748405504596428952168554d0];  end;
if firstCall,   wg(15)=[0.102852652893558840341285636705415d0];  end;
%
if firstCall,   xgk(1)=[0.999484410050490637571325895705811d0];  end;
if firstCall,   xgk(2)=[0.996893484074649540271630050918695d0];  end;
if firstCall,   xgk(3)=[0.991630996870404594858628366109486d0];  end;
if firstCall,   xgk(4)=[0.983668123279747209970032581605663d0];  end;
if firstCall,   xgk(5)=[0.973116322501126268374693868423707d0];  end;
if firstCall,   xgk(6)=[0.960021864968307512216871025581798d0];  end;
if firstCall,   xgk(7)=[0.944374444748559979415831324037439d0];  end;
if firstCall,   xgk(8)=[0.926200047429274325879324277080474d0];  end;
if firstCall,   xgk(9)=[0.905573307699907798546522558925958d0];  end;
if firstCall,   xgk(10)=[0.882560535792052681543116462530226d0];  end;
if firstCall,   xgk(11)=[0.857205233546061098958658510658944d0];  end;
if firstCall,   xgk(12)=[0.829565762382768397442898119732502d0];  end;
if firstCall,   xgk(13)=[0.799727835821839083013668942322683d0];  end;
if firstCall,   xgk(14)=[0.767777432104826194917977340974503d0];  end;
if firstCall,   xgk(15)=[0.733790062453226804726171131369528d0];  end;
if firstCall,   xgk(16)=[0.697850494793315796932292388026640d0];  end;
if firstCall,   xgk(17)=[0.660061064126626961370053668149271d0];  end;
if firstCall,   xgk(18)=[0.620526182989242861140477556431189d0];  end;
if firstCall,   xgk(19)=[0.579345235826361691756024932172540d0];  end;
if firstCall,   xgk(20)=[0.536624148142019899264169793311073d0];  end;
if firstCall,   xgk(21)=[0.492480467861778574993693061207709d0];  end;
if firstCall,   xgk(22)=[0.447033769538089176780609900322854d0];  end;
if firstCall,   xgk(23)=[0.400401254830394392535476211542661d0];  end;
if firstCall,   xgk(24)=[0.352704725530878113471037207089374d0];  end;
if firstCall,   xgk(25)=[0.304073202273625077372677107199257d0];  end;
if firstCall,   xgk(26)=[0.254636926167889846439805129817805d0];  end;
if firstCall,   xgk(27)=[0.204525116682309891438957671002025d0];  end;
if firstCall,   xgk(28)=[0.153869913608583546963794672743256d0];  end;
if firstCall,   xgk(29)=[0.102806937966737030147096751318001d0];  end;
if firstCall,   xgk(30)=[0.051471842555317695833025213166723d0];  end;
if firstCall,   xgk(31)=[0.000000000000000000000000000000000d0];  end;
%
if firstCall,   wgk(1)=[0.001389013698677007624551591226760d0];  end;
if firstCall,   wgk(2)=[0.003890461127099884051267201844516d0];  end;
if firstCall,   wgk(3)=[0.006630703915931292173319826369750d0];  end;
if firstCall,   wgk(4)=[0.009273279659517763428441146892024d0];  end;
if firstCall,   wgk(5)=[0.011823015253496341742232898853251d0];  end;
if firstCall,   wgk(6)=[0.014369729507045804812451432443580d0];  end;
if firstCall,   wgk(7)=[0.016920889189053272627572289420322d0];  end;
if firstCall,   wgk(8)=[0.019414141193942381173408951050128d0];  end;
if firstCall,   wgk(9)=[0.021828035821609192297167485738339d0];  end;
if firstCall,   wgk(10)=[0.024191162078080601365686370725232d0];  end;
if firstCall,   wgk(11)=[0.026509954882333101610601709335075d0];  end;
if firstCall,   wgk(12)=[0.028754048765041292843978785354334d0];  end;
if firstCall,   wgk(13)=[0.030907257562387762472884252943092d0];  end;
if firstCall,   wgk(14)=[0.032981447057483726031814191016854d0];  end;
if firstCall,   wgk(15)=[0.034979338028060024137499670731468d0];  end;
if firstCall,   wgk(16)=[0.036882364651821229223911065617136d0];  end;
if firstCall,   wgk(17)=[0.038678945624727592950348651532281d0];  end;
if firstCall,   wgk(18)=[0.040374538951535959111995279752468d0];  end;
if firstCall,   wgk(19)=[0.041969810215164246147147541285970d0];  end;
if firstCall,   wgk(20)=[0.043452539701356069316831728117073d0];  end;
if firstCall,   wgk(21)=[0.044814800133162663192355551616723d0];  end;
if firstCall,   wgk(22)=[0.046059238271006988116271735559374d0];  end;
if firstCall,   wgk(23)=[0.047185546569299153945261478181099d0];  end;
if firstCall,   wgk(24)=[0.048185861757087129140779492298305d0];  end;
if firstCall,   wgk(25)=[0.049055434555029778887528165367238d0];  end;
if firstCall,   wgk(26)=[0.049795683427074206357811569379942d0];  end;
if firstCall,   wgk(27)=[0.050405921402782346840893085653585d0];  end;
if firstCall,   wgk(28)=[0.050881795898749606492297473049805d0];  end;
if firstCall,   wgk(29)=[0.051221547849258772170656282604944d0];  end;
if firstCall,   wgk(30)=[0.051426128537459025933862879215781d0];  end;
if firstCall,   wgk(31)=[0.051494729429451567558340433647099d0];  end;
firstCall=0;
%
%           LIST OF MAJOR VARIABLES
%           -----------------------
%
%           CENTR  - MID POINT OF THE INTERVAL
%           HLGTH  - HALF-LENGTH OF THE INTERVAL
%           DABSC  - ABSCISSA
%           FVAL*  - FUNCTION VALUE
%           RESG   - RESULT OF THE 30-POINT GAUSS RULE
%           RESK   - RESULT OF THE 61-POINT KRONROD RULE
%           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F
%                    OVER (A,B), I.E. TO I/(B-A)
%
%           MACHINE DEPENDENT CONSTANTS
%           ---------------------------
%
%           EPMACH IS THE LARGEST RELATIVE SPACING.
%           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
%
%***FIRST EXECUTABLE STATEMENT  DQK61
[epmach ]=d1mach(4);
[uflow ]=d1mach(1);
%
centr = 0.5d+00.*(b+a);
hlgth = 0.5d+00.*(b-a);
dhlgth = abs(hlgth);
%
%           COMPUTE THE 61-POINT KRONROD APPROXIMATION TO THE
%           INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
%
resg = 0.0d+00;
fc = f(centr);
resk = wgk(31).*fc;
resabs = abs(resk);
for j = 1 : 15;
jtw = fix(j.*2);
dabsc = hlgth.*xgk(jtw);
fval1 = f(centr-dabsc);
fval2 = f(centr+dabsc);
fv1(jtw) = fval1;
fv2(jtw) = fval2;
fsum = fval1 + fval2;
resg = resg + wg(j).*fsum;
resk = resk + wgk(jtw).*fsum;
resabs = resabs + wgk(jtw).*(abs(fval1)+abs(fval2));
end; j = fix(15+1);
for j = 1 : 15;
jtwm1 = fix(j.*2 - 1);
dabsc = hlgth.*xgk(jtwm1);
fval1 = f(centr-dabsc);
fval2 = f(centr+dabsc);
fv1(jtwm1) = fval1;
fv2(jtwm1) = fval2;
fsum = fval1 + fval2;
resk = resk + wgk(jtwm1).*fsum;
resabs = resabs + wgk(jtwm1).*(abs(fval1)+abs(fval2));
end; j = fix(15+1);
reskh = resk.*0.5d+00;
resasc = wgk(31).*abs(fc-reskh);
for j = 1 : 30;
resasc = resasc + wgk(j).*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh));
end; j = fix(30+1);
result = resk.*hlgth;
resabs = resabs.*dhlgth;
resasc = resasc.*dhlgth;
abserr = abs((resk-resg).*hlgth);
if( resasc~=0.0d+00 && abserr~=0.0d+00 )
abserr = resasc.*min(0.1d+01,(0.2d+03.*abserr./resasc).^1.5d+00);
end;
if( resabs>uflow./(0.5d+02.*epmach) )
abserr = max((epmach.*0.5d+02).*resabs,abserr);
end;
end
%DECK DQMOMO
