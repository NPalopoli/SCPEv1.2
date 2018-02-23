/*calculate eigensystem - functions  were taken from MOLPHY  Jun Adachi (C) 1996*/
/*ProtML 2.3b3(July 1 1996) Maximum Likelihood Inference of Protein Phylogeny
Copyright (C) 1992-1996 J. Adachi & M. Hasegawa. All rights reserved.*/


#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "scpe.h"

#define CUAA 8000
#define SQ 400

double Q[AA][AA];
int vector[AA];
double Eval[AA];
double evali[AA];
double Evec[AA][AA];
double invEvec[AA][AA], dummy[AA][AA];
int err;
double Cijk[CUAA];
double lambda,qii[AA];



int diag(double *qmatrix,double *pmatrix,double time)
{
	
	int i,j,k;
        float sum;

	for (i=0; i<AA; i++) for (j=0; j<AA; j++) Q[i][j]=*(qmatrix+(i*AA+j));
	for (i=0; i<AA; i++) for (j=0; j<AA; j++) *(pmatrix+(i*AA+j))=0.;
	


		
	elmhes(Q, vector, AA);			   
	eltran(Q, Evec, vector, AA);	
	hqr2(AA, 1, AA, &err, Q, Evec, Eval, evali);
	for(i=0;i<AA;i++) for(j=0;j<AA;j++) dummy[i][j]=Evec[i][j];
   /* printEigen();*/
	luinverse(dummy,invEvec,AA);
		
	for (i=0; i<AA; i++){
		for (j=0; j<AA; j++){
			for (k=0; k<AA; k++)
					Cijk[i*AA*AA+j*AA+k] = Evec[i][k]*invEvec[k][j];
			}
	}
		
	SetMatrix(pmatrix,time);
      /*print equilibrium frequences*/

/*	printf("Frecuencias\n");
	for (i=0; i<AA; i++) {
                for (j=0; j<AA; j++) printf("%f ",*(pmatrix+(i*AA+j)));
                printf("\n");
        }*/
    


	return(0);
}

void printEigen()
{
int i, j;

	fprintf(stderr, "\n\nEigensystem:\nVectors:\n");
	for(i=0;i<AA;i++){
		for(j=0;j<AA;j++)
			fprintf(stderr, "%lf ", Evec[i][j]);
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "\nValues: ");
	for(j=0;j<AA;j++)
		fprintf(stderr, "%lf ", Eval[j]);
	fprintf(stderr, "\n");
}

void elmhes(dmattpmty a, ivectpmty ordr, int n)
{
	int m, j, i;
	double y, x;

	for (i = 0; i < n; i++)
		ordr[i] = 0;
	for (m = 2; m < n; m++) {
		x = 0.0;
		i = m;
		for (j = m; j <= n; j++) {
			if (fabs(a[j - 1][m - 2]) > fabs(x)) {
				x = a[j - 1][m - 2];
				i = j;
			}
		}
		ordr[m - 1] = i;      /* vector */
		if (i != m) {
			for (j = m - 2; j < n; j++) {
				y = a[i - 1][j];
				a[i - 1][j] = a[m - 1][j];
				a[m - 1][j] = y;
			}
			for (j = 0; j < n; j++) {
				y = a[j][i - 1];
				a[j][i - 1] = a[j][m - 1];
				a[j][m - 1] = y;
			}
		}
		if (x != 0.0) {
			for (i = m; i < n; i++) {
				y = a[i][m - 2];
				if (y != 0.0) {
					y /= x;
					a[i][m - 2] = y;
					for (j = m - 1; j < n; j++)
						a[i][j] -= y * a[m - 1][j];
					for (j = 0; j < n; j++)
						a[j][m - 1] += y * a[j][i];
				}
			}
		}
	}
} /*_ elmhes */

void eltran(dmattpmty a, dmattpmty zz, ivectpmty ordr, int n)
{
	int i, j, m;

	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			zz[i][j] = 0.0;
			zz[j][i] = 0.0;
		}
		zz[i][i] = 1.0;
	}
	if (n <= 2)
		return;
	for (m = n - 1; m >= 2; m--) {
		for (i = m; i < n; i++)
			zz[i][m - 1] = a[i][m - 2];
		i = ordr[m - 1];
		if (i != m) {
			for (j = m - 1; j < n; j++) {
				zz[m - 1][j] = zz[i - 1][j];
				zz[i - 1][j] = 0.0;
			}
			zz[i - 1][m - 1] = 1.0;
		}
	}
} /*_ eltran */

void hqr2(int n, int low, int hgh, int *err, dmattpmty h, dmattpmty zz, dvectpmty wr, dvectpmty wi)
{
	int i, j, k, l, m, en, na, itn, its;
	double p, q, r, s, t, w, x, y, ra, sa, vi, vr, z, norm, tst1, tst2;
	boolean notlas;

	*err = 0;
	norm = 0.0;
	k = 1;
	/* store roots isolated by balanc and compute matrix norm */
	for (i = 0; i < n; i++) {
		for (j = k - 1; j < n; j++)
			norm += fabs(h[i][j]);
		k = i + 1;
		if (i + 1 < low || i + 1 > hgh) {
			wr[i] = h[i][i];
			wi[i] = 0.0;
		}
	}
	en = hgh;
	t = 0.0;
	itn = n * 30;

	while (en >= low) {	       /* search for next eigenvalues */
		its = 0;
		na = en - 1;

		while (en >= 1) {      /* infinietr loop */
			/* look for single small sub-diagonal element */
			for (l = en; l > low; l--) {
				s = fabs(h[l - 2][l - 2]) + fabs(h[l - 1][l - 1]);
				if (s == 0.0)
					s = norm;
				tst1 = s;
				tst2 = tst1 + fabs(h[l - 1][l - 2]);
				if (tst2 == tst1)
					goto L100;
			}
			l = low;
	L100:
			x = h[en - 1][en - 1];	/* form shift */
			if (l == en || l == na)
				break;
			if (itn == 0) { /* all eigenvalues have not converged */
				*err = en;
				goto Lerror;
			}
			y = h[na - 1][na - 1];
			w = h[en - 1][na - 1] * h[na - 1][en - 1];
			/* form exceptional shift */
			if (its == 10 || its == 20) {
				t += x;
				for (i = low - 1; i < en; i++)
					h[i][i] -= x;
				s = fabs(h[en - 1][na - 1]) + fabs(h[na - 1][en - 3]);
				x = 0.75 * s;
				y = x;
				w = -0.4375 * s * s;
			}
			its++;
			itn--;
			/* look for two consecutive small sub-diagonal elements */
			for (m = en - 2; m >= l; m--) {
				z = h[m - 1][m - 1];
				r = x - z;
				s = y - z;
				p = (r * s - w) / h[m][m - 1] + h[m - 1][m];
				q = h[m][m] - z - r - s;
				r = h[m + 1][m];
				s = fabs(p) + fabs(q) + fabs(r);
				p /= s;
				q /= s;
				r /= s;
				if (m == l)
					break;
				tst1 = fabs(p) * (fabs(h[m-2][m-2]) + fabs(z) + fabs(h[m][m]));
				tst2 = tst1 + fabs(h[m - 1][m - 2]) * (fabs(q) + fabs(r));
				if (tst2 == tst1)
					break;
			}

			for (i = m + 2; i <= en; i++) {
				h[i - 1][i - 3] = 0.0;
				if (i != m + 2)
					h[i - 1][i - 4] = 0.0;
			}

			/* double qr step involving rows l to en and columns m to en */
			for (k = m; k <= na; k++) {
				notlas = (k != na);
				if (k != m) {
					p = h[k - 1][k - 2];
					q = h[k][k - 2];
					r = 0.0;
					if (notlas)
						r = h[k + 1][k - 2];
					x = fabs(p) + fabs(q) + fabs(r);
					if (x != 0.0) {
						p /= x;
						q /= x;
						r /= x;
					}
				}
				if (x != 0.0) {
					if (p < 0.0) /* sign */
						s = - sqrt(p * p + q * q + r * r);
					else
						s = sqrt(p * p + q * q + r * r);
					if (k != m)
						h[k - 1][k - 2] = -s * x;
					else {
						if (l != m)
							h[k - 1][k - 2] = -h[k - 1][k - 2];
					}
					p += s;
					x = p / s;
					y = q / s;
					z = r / s;
					q /= p;
					r /= p;
					if (!notlas) {
						for (j = k - 1; j < n; j++) {	/* row modification */
							p = h[k - 1][j] + q * h[k][j];
							h[k - 1][j] -= p * x;
							h[k][j] -= p * y;
						}
						j = (en < (k + 3)) ? en : (k + 3); /* min */
						for (i = 0; i < j; i++) {	/* column modification */
							p = x * h[i][k - 1] + y * h[i][k];
							h[i][k - 1] -= p;
							h[i][k] -= p * q;
						}
						/* accumulate transformations */
						for (i = low - 1; i < hgh; i++) {
							p = x * zz[i][k - 1] + y * zz[i][k];
							zz[i][k - 1] -= p;
							zz[i][k] -= p * q;
						}
					} else {
						for (j = k - 1; j < n; j++) {	/* row modification */
							p = h[k - 1][j] + q * h[k][j] + r * h[k + 1][j];
							h[k - 1][j] -= p * x;
							h[k][j] -= p * y;
							h[k + 1][j] -= p * z;
						}
						j = (en < (k + 3)) ? en : (k + 3); /* min */
						for (i = 0; i < j; i++) {	/* column modification */
							p = x * h[i][k - 1] + y * h[i][k] + z * h[i][k + 1];
							h[i][k - 1] -= p;
							h[i][k] -= p * q;
							h[i][k + 1] -= p * r;
						}
						/* accumulate transformations */
						for (i = low - 1; i < hgh; i++) {
							p = x * zz[i][k-1] + y * zz[i][k] + z * zz[i][k+1];
							zz[i][k - 1] -= p;
							zz[i][k] -= p * q;
							zz[i][k + 1] -= p * r;
						}
					}
				}
			}	       /* for k */
		}		       /* while infinite loop */


		if (l == en) {	       /* one root found */
			h[en - 1][en - 1] = x + t;
			wr[en - 1] = h[en - 1][en - 1];
			wi[en - 1] = 0.0;
			en = na;
			continue;
		}
		y = h[na - 1][na - 1];
		w = h[en - 1][na - 1] * h[na - 1][en - 1];
		p = (y - x) / 2.0;
		q = p * p + w;
		z = sqrt(fabs(q));
		h[en - 1][en - 1] = x + t;
		x = h[en - 1][en - 1];
		h[na - 1][na - 1] = y + t;
		if (q >= 0.0) {	       /* real pair */
			if (p < 0.0) /* sign */
				z = p - fabs(z);
			else
				z = p + fabs(z);
			wr[na - 1] = x + z;
			wr[en - 1] = wr[na - 1];
			if (z != 0.0)
				wr[en - 1] = x - w / z;
			wi[na - 1] = 0.0;
			wi[en - 1] = 0.0;
			x = h[en - 1][na - 1];
			s = fabs(x) + fabs(z);
			p = x / s;
			q = z / s;
			r = sqrt(p * p + q * q);
			p /= r;
			q /= r;
			for (j = na - 1; j < n; j++) {	/* row modification */
				z = h[na - 1][j];
				h[na - 1][j] = q * z + p * h[en - 1][j];
				h[en - 1][j] = q * h[en - 1][j] - p * z;
			}
			for (i = 0; i < en; i++) {	/* column modification */
				z = h[i][na - 1];
				h[i][na - 1] = q * z + p * h[i][en - 1];
				h[i][en - 1] = q * h[i][en - 1] - p * z;
			}
			/* accumulate transformations */
			for (i = low - 1; i < hgh; i++) {
				z = zz[i][na - 1];
				zz[i][na - 1] = q * z + p * zz[i][en - 1];
				zz[i][en - 1] = q * zz[i][en - 1] - p * z;
			}
		} else {	       /* complex pair */
			wr[na - 1] = x + p;
			wr[en - 1] = x + p;
			wi[na - 1] = z;
			wi[en - 1] = -z;
		}
		en -= 2;
	}			       /* while en >= low */

	/* backsubstitute to find vectors of upper triangular form */
	if (norm != 0.0) {
		for (en = n; en >= 1; en--) {
			p = wr[en - 1];
			q = wi[en - 1];
			na = en - 1;
			if (q == 0.0) {/* real vector */
				m = en;
				h[en - 1][en - 1] = 1.0;
				if (na != 0) {
					for (i = en - 2; i >= 0; i--) {
						w = h[i][i] - p;
						r = 0.0;
						for (j = m - 1; j < en; j++)
							r += h[i][j] * h[j][en - 1];
						if (wi[i] < 0.0) {
							z = w;
							s = r;
						} else {
							m = i + 1;
							if (wi[i] == 0.0) {
								t = w;
								if (t == 0.0) {
									tst1 = norm;
									t = tst1;
									do {
										t = 0.01 * t;
										tst2 = norm + t;
									} while (tst2 > tst1);
								}
								h[i][en - 1] = -(r / t);
							} else {	/* solve real equations */
								x = h[i][i + 1];
								y = h[i + 1][i];
								q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
								t = (x * s - z * r) / q;
								h[i][en - 1] = t;
								if (fabs(x) > fabs(z))
									h[i + 1][en - 1] = (-r - w * t) / x;
								else
									h[i + 1][en - 1] = (-s - y * t) / z;
							}
							/* overflow control */
							t = fabs(h[i][en - 1]);
							if (t != 0.0) {
								tst1 = t;
								tst2 = tst1 + 1.0 / tst1;
								if (tst2 <= tst1) {
									for (j = i; j < en; j++)
										h[j][en - 1] /= t;
								}
							}
						}
					}
				}
			} else if (q > 0.0) {
				m = na;
				if (fabs(h[en - 1][na - 1]) > fabs(h[na - 1][en - 1])) {
					h[na - 1][na - 1] = q / h[en - 1][na - 1];
					h[na - 1][en - 1] = (p - h[en-1][en-1]) / h[en-1][na-1];
				} else
					cdiv(0.0, -h[na-1][en-1], h[na-1][na-1] - p, q, &h[na-1][na-1], &h[na-1][en-1]);
				h[en - 1][na - 1] = 0.0;
				h[en - 1][en - 1] = 1.0;
				if (en != 2) {
					for (i = en - 3; i >= 0; i--) {
						w = h[i][i] - p;
						ra = 0.0;
						sa = 0.0;
						for (j = m - 1; j < en; j++) {
							ra += h[i][j] * h[j][na - 1];
							sa += h[i][j] * h[j][en - 1];
						}
						if (wi[i] < 0.0) {
							z = w;
							r = ra;
							s = sa;
						} else {
							m = i + 1;
							if (wi[i] == 0.0)
								cdiv(-ra, -sa, w, q, &h[i][na-1], &h[i][en-1]);
							else {	/* solve complex equations */
								x = h[i][i + 1];
								y = h[i + 1][i];
								vr = (wr[i]-p) * (wr[i]-p) + wi[i]*wi[i] - q*q;
								vi = (wr[i] - p) * 2.0 * q;
								if (vr == 0.0 && vi == 0.0) {
									tst1 = norm * (fabs(w) + fabs(q) + fabs(x)
										+ fabs(y) + fabs(z));
									vr = tst1;
									do {
										vr = 0.01 * vr;
										tst2 = tst1 + vr;
									} while (tst2 > tst1);
								}
								cdiv(x * r - z * ra + q * sa,
									x * s - z * sa - q * ra, vr, vi,
								     &h[i][na - 1], &h[i][en - 1]);
								if (fabs(x) > fabs(z) + fabs(q)) {
									h[i + 1][na - 1] = (q * h[i][en - 1]
										- w * h[i][na - 1] - ra) / x;
									h[i + 1][en - 1] = (-sa - w * h[i][en - 1]
										- q * h[i][na - 1]) / x;
								} else
									cdiv(-r - y * h[i][na - 1],
										-s - y * h[i][en - 1], z, q,
									     &h[i + 1][na - 1], &h[i + 1][en - 1]);
							}
							/* overflow control */
							t = (fabs(h[i][na-1]) > fabs(h[i][en-1])) ?
								 fabs(h[i][na-1]) : fabs(h[i][en-1]); /* max */
							if (t != 0.0) {
								tst1 = t;
								tst2 = tst1 + 1.0 / tst1;
								if (tst2 <= tst1) {
									for (j = i; j < en; j++) {
										h[j][na - 1] /= t;
										h[j][en - 1] /= t;
									}
								}
							}
						}
					}
				}
			}
		}		       /* for nn */

		/* end back substitution. vectors of isolated roots */
		for (i = 0; i < n; i++) {
			if (i + 1 < low || i + 1 > hgh) {
				for (j = i; j < n; j++)
					zz[i][j] = h[i][j];
			}
		}
		/* multiply by transformation matrix to give vectors of
		 * original full matrix. */
		for (j = n - 1; j >= low - 1; j--) {
			m = ((j + 1) < hgh) ? (j + 1) : hgh; /* min */
			for (i = low - 1; i < hgh; i++) {
				z = 0.0;
				for (k = low - 1; k < m; k++)
					z += zz[i][k] * h[k][j];
				zz[i][j] = z;
			}
		}
	}
	return;

Lerror: /* two roots found and complex vector */
	fprintf(stderr, "ERROR in hqr2. two roots found or complex vector");
	exit(1);

} /*_ hqr2 */

static void cdiv(double ar,double ai,double br,double bi,double *cr,double *ci)
{
	double s, ars, ais, brs, bis;

	s = fabs(br) + fabs(bi);
	ars = ar / s;
	ais = ai / s;
	brs = br / s;
	bis = bi / s;
	s = brs * brs + bis * bis;
	*cr = (ars * brs + ais * bis) / s;
	*ci = (ais * brs - ars * bis) / s;
} /*_ cdiv */
	
void luinverse(dmattpmty omtrx, dmattpmty imtrx, int size)
{
	/* INVERSION OF MATRIX ON LU DECOMPOSITION */
    double eps = 1.0e-20; /* ! */
	int i, j, k, l, maxi, idx, ix, jx;
	double sum, tmp, maxb, aw;
	ivectpmty index;
	double *wk;

	wk = (double *) malloc((unsigned)size * sizeof(double));
	aw = 1.0;
	for (i = 0; i < size; i++) {
		maxb = 0.0;
		for (j = 0; j < size; j++) {
			if (fabs(omtrx[i][j]) > maxb)
				maxb = fabs(omtrx[i][j]);
		}
		if (maxb == 0.0) {
			fprintf(stderr, "luinverse: singular matrix\n");
			exit(1);
		}
		wk[i] = 1.0 / maxb;
	}
	for (j = 0; j < size; j++) {
		for (i = 0; i < j; i++) {
			sum = omtrx[i][j];
			for (k = 0; k < i; k++)
				sum -= omtrx[i][k] * omtrx[k][j];
			omtrx[i][j] = sum;
		}
		maxb = 0.0;
		for (i = j; i < size; i++) {
			sum = omtrx[i][j];
			for (k = 0; k < j; k++)
				sum -= omtrx[i][k] * omtrx[k][j];
			omtrx[i][j] = sum;
			tmp = wk[i] * fabs(sum);
			if (tmp >= maxb) {
				maxb = tmp;
				maxi = i;
			}
		}
		if (j != maxi) {
			for (k = 0; k < size; k++) {
				tmp = omtrx[maxi][k];
				omtrx[maxi][k] = omtrx[j][k];
				omtrx[j][k] = tmp;
			}
			aw = -aw;
			wk[maxi] = wk[j];
		}
		index[j] = maxi;
		if (omtrx[j][j] == 0.0)
			omtrx[j][j] = eps;
		if (j != size - 1) {
			tmp = 1.0 / omtrx[j][j];
			for (i = j + 1; i < size; i++)
				omtrx[i][j] *= tmp;
		}
	}
	for (jx = 0; jx < size; jx++) {
		for (ix = 0; ix < size; ix++)
			wk[ix] = 0.0;
		wk[jx] = 1.0;
		l = -1;
		for (i = 0; i < size; i++) {
			idx = index[i];
			sum = wk[idx];
			wk[idx] = wk[i];
			if (l != -1) {
				for (j = l; j < i; j++)
					sum -= omtrx[i][j] * wk[j];
			} else if (sum != 0.0)
				l = i;
			wk[i] = sum;
		}
		for (i = size - 1; i >= 0; i--) {
			sum = wk[i];
			for (j = i + 1; j < size; j++)
				sum -= omtrx[i][j] * wk[j];
			wk[i] = sum / omtrx[i][i];
		}
		for (ix = 0; ix < size; ix++)
			imtrx[ix][jx] = wk[ix];
	}
	/*free((char *)wk);gustavo parisi*/
	wk = NULL;
} /*_ luinverse */


void SetMatrix(double *matrix, double len)
{	
	int i,j,k;
	double expt[AA];
	double *P, *Q;
	
/* P(t)ij = SUM Cijk * exp{Root*t}
*/
		P=matrix;
		if (len<1e-6) { 
			for (i=0; i<AA; i++) 
				for (j=0; j<AA; j++) {
					if (i==j)
						*P=1.0;
					else 	
						*P=0.0;
					P++;
				}
			return; 
		}
		
		for (k=0; k<AA; k++) expt[k]=exp(len*Eval[k]);
		for (i=0; i<AA; i++) 
			for (j=0; j<AA; j++) {
				(*P)=0.0;
				for (k=0; k<AA; k++)
					(*P)+=Cijk[i*AA*AA+j*AA+k]*expt[k];
				P++;
			}



/*	P=Q=matrix; make matrix cumulative      G. Parisi 2002
	for(i=0;i<AA;i++){
		for(j=1;j<AA;j++){
			P++;
			*P+=(*Q);
			Q++;
		}
		P++;
		Q++;
	}*/
}


