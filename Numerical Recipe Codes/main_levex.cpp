#include <iostream>
#include <cmath>
#include <iomanip>
#include "series.h"
#include "bessel.h"
#include "levex.h"

using namespace std;
int main_levex(void)
{
const Doub PI=3.141592653589793;
int nterm=12;
doub beta=1.0,a=0.0,b=0.0,sum=0.0;
Levin series(100,0.0);
cout << setw(5) << "N" << setw(19) << "Sum (direct)" << setw(21)
<< "Sum (Levin)" << endl;
for (Int n=0; n<=nterm; n++) {
b+=PI;
doub s=qromb(func,a,b,1.e-8);
a=b;
sum+=s;
doub omega=(beta+n)*s; Use u transformation.
doub ans=series.next(sum,omega,beta);
cout << setw(5) << n << fixed << setprecision(14) << setw(21)
<< sum << setw(21) << ans << endl;
}
return 0;
}
