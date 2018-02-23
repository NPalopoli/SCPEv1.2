/******************************************************************************
 random.h
 
 
   U(0,1): AS 183: Appl. Stat. 31:188-190 
   Wichmann BA & Hill ID.  1982.  An efficient and portable
   pseudo-random number generator.  Appl. Stat. 31:188-190

   x, y, z are any numbers in the range 1-30000.  Integer operation up
   to 30323 required.
 	
******************************************************************************/

#ifndef _RANDOM_H_
#define _RANDOM_H_

void SetSeed (int seed);
double rndu (void);
double rndgamma (double s);

#endif /* _RANDOM_H_ */
