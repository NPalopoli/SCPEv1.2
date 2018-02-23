/******************************************************************************
 
   The following code comes from tools.c in Yang's PAML package

   U(0,1): AS 183: Appl. Stat. 31:188-190 
   Wichmann BA & Hill ID.  1982.  An efficient and portable
   pseudo-random number generator.  Appl. Stat. 31:188-190

   x, y, z are any numbers in the range 1-30000.  Integer operation up
   to 30323 required.
 	
******************************************************************************/


#include <stdio.h>
#include <math.h>
#include "random.h"


static int z_rndu=137;
static unsigned w_rndu=13757;

void SetSeed (int seed)
{
   z_rndu = 170*(seed%178) + 137;
   w_rndu=seed;
}


#ifdef FAST_RANDOM_NUMBER

double rndu (void)
{
   w_rndu *= 127773;
   return ldexp((double)w_rndu, -32);
}

#else 

double rndu (void)
{
   static int x_rndu=11, y_rndu=23;
   double r;

   x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
   y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
   z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
   if (x_rndu<0) x_rndu+=30269;
   if (y_rndu<0) y_rndu+=30307;
   if (z_rndu<0) z_rndu+=30323;
   r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
   return (r-(int)r);
}

#endif
double rndgamma1 (double s);
double rndgamma2 (double s);

double rndgamma (double s)
{
	double	r=0.0;
	
	if (s <= 0.0)      
		puts ("Error gamma..");
	else if (s < 1.0)  
		r = rndgamma1 (s);
	else if (s > 1.0)  
		r = rndgamma2 (s);
	else           
		r =- log(rndu());
	return (r);
}


double rndgamma1 (double s)
{

	double			r, x=0.0, small=1e-37, w;
	static double	a, p, uf, ss=10.0, d;
	
	if (s!=ss) 
		{
		a  = 1.0-s;
		p  = a/(a+s*exp(-a));
		uf = p*pow(small/a,s);
		d  = a*log(a);
		ss = s;
		}
	for (;;) 
		{
		r = rndu();
		if (r > p)        
			x = a-log((1.0-r)/(1.0-p)), w=a*log(x)-d;
		else if (r>uf)  
			x = a*pow(r/p,1/s), w=x;
		else            
			return (0.0);
		r = rndu();
		if (1.0-r <= w && r > 0.0)
			if (r*(w+1.0) >= 1.0 || -log(r) <= w)  
				continue;
		break;
		}
	return (x);
}


double rndgamma2 (double s)
{

	double			r ,d, f, g, x;
	static double	b, h, ss=0;
	
	if (s!=ss) 
		{
		b  = s-1.0;
		h  = sqrt(3.0*s-0.75);
		ss = s;
		}
	for (;;) 
		{
		r = rndu();
		g = r-r*r;
		f = (r-0.5)*h/sqrt(g);
		x = b+f;
		if (x <= 0.0) 
			continue;
		r = rndu();
		d = 64*r*r*g*g*g;
		if (d*x < x-2.0*f*f || log(d) < 2*(b*log(x/b)-f))  
			break;
		}
	return (x);
}



