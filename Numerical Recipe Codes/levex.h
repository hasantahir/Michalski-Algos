Doub func(const Doub x)
{
	if (x == 0.0)
		return 0.0;
	else {
		Bessel bess;
		return x*bess.jnu(0.0,x)/(1.0+x*x);
	}
}
