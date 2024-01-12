/*
	cuhre-f.c
		Fortran interface for Cuhre
		this file is part of Cuhre
		last modified 9 Feb 05 th
*/


#include "decl.h"

/*
#ifdef HAVE_UNDERSCORE
#define cuhre cuhre_
#endif
*/

Extern void Cuhre(ccount ndim, ccount ncomp, Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, ccount mineval, ccount maxeval,
  ccount key,
  count *pnregions, count *pneval, int *pfail,
  real *integral, real *error, real *prob);


Extern void cuhre_(ccount *pndim, ccount *pncomp, Integrand integrand,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, ccount *pmineval, ccount *pmaxeval,
  ccount *pkey,
  count *pnregions, count *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  Cuhre(*pndim, *pncomp, integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pkey,
    pnregions, pneval, pfail,
    integral, error, prob);
}

