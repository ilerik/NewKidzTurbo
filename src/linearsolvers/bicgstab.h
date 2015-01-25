#include "basetypes.h"

// -*- C++ -*-

#ifndef BICGSTAB_BLAS_H
#define BICGSTAB_BLAS_H

// ============================================================================
//
//  BICGstab
//
//  siehe
//  Ashby, Manteuffel, Saylor
//     A taxononmy for conjugate gradient methods
//     SIAM J Numer Anal 27, 1542-1568 (1990)
//
//  oder:
//  Willy D"orfler:
//     Orthogonale Fehlermethoden
//
//                                                 ----------------------------
//                                                 Christian Badura, Mai 1998
//
// ============================================================================



template< class MATRIX > inline
int
bicgstab( unsigned N, const MATRIX &A,
	  const double *b, double *x, double eps );


template< class MATRIX > inline
int
bicgstab( unsigned N, const MATRIX &A,
	  const double *b, double *x, double eps, bool detailed );



template< class MATRIX, class PC_MATRIX > inline
int
bicgstab( unsigned N, const MATRIX &A, const PC_MATRIX &C,
	  const double *b, double *x, double eps );


template< class MATRIX, class PC_MATRIX > inline
int
bicgstab( unsigned N, const MATRIX &A, const PC_MATRIX &C,
	  const double *b, double *x, double eps, bool detailed );


// ============================================================================


#include <cassert>
//#include "cblas.h"
#include "mkl_blas.h"

// ============================================================================


template< class MATRIX > inline
int
bicgstab( unsigned NN, MATRIX &A,
	  const double *b, double *x, double eps, bool detailed ) {
	int N = NN;	
  double *rT  = new double[N];
  double *d   = new double[N];
  double *h   = new double[N];
  double *u   = new double[N];
  double *Ad  = new double[N];
  double *t   = new double[N];
  double *s   = h;
  double rTh, rTAd, rTr, alpha, beta, omega, st, tt;
  int its=0;
  int one=1;
  double err=eps*eps*ddot(&N,b,&one,b,&one);
  // f"ur's Abbruchkriterium (*)  -- r enth"alt immer das Residuum r=Ax-b
  double *r = new double[N];

  mult(A,x,r);
  double minusone = -1.0;
  daxpy(&N, &minusone,b,&one,r,&one);
  dcopy(&N,r,&one,d,&one);
  dcopy(&N,d,&one,h,&one);
  dcopy(&N,h,&one,rT,&one);
  assert( ddot(&N,rT,&one,rT,&one)>1e-40 );
  rTh=ddot(&N,rT,&one,h,&one);
  rTr=ddot(&N,r,&one,r,&one);
  while ( rTr>err ) {
    mult(A,d,Ad);
    rTAd=ddot(&N,rT,&one,Ad,&one);
    assert( fabs(rTAd)>1e-40 );
    alpha=rTh/rTAd;
	double minusalpha = -alpha;
    daxpy(&N,&minusalpha,Ad,&one,r,&one);
    dcopy(&N,h,&one,s,&one);
    daxpy(&N,&minusalpha,Ad,&one,s,&one);
    mult(A,s,t);
	double oned = 1.0;
    daxpy(&N,&oned,t,&one,u,&one);
    dscal(&N,&alpha,u,&one);
    st=ddot(&N,s,&one,t,&one);
    tt=ddot(&N,t,&one,t,&one);
    if ( fabs(st)<1e-40 || fabs(tt)<1e-40 )
      omega = 0.;
    else
      omega = st/tt;
	double minusomega = -omega;
	minusalpha = -alpha;
    daxpy(&N, &minusomega,t,&one,r,&one);
    daxpy(&N, &minusalpha,d,&one,x,&one);
    daxpy(&N, &minusomega,s,&one,x,&one);
    dcopy(&N,s,&one,h,&one);
    daxpy(&N, &minusomega,t,&one,h,&one);
    beta=(alpha/omega)/rTh; rTh=ddot(&N,rT,&one,h,&one); beta*=rTh;
    dscal(&N,&beta,d,&one);
    daxpy(&N,&oned,h,&one,d,&one);
	double minusbetaomega = minusomega * beta;
    daxpy(&N, &minusbetaomega,Ad,&one,d,&one);
    rTr=ddot(&N,r,&one,r,&one);
    if ( detailed )
      std::cout<<"bicgstab "<<its<<"\t"<<sqrt(rTr)<<std::endl;
    ++its;
  }
  delete[] r;
  delete[] rT;
  delete[] d;
  delete[] h;
  delete[] u;
  delete[] Ad;
  delete[] t;
  return its;
}


template< class MATRIX > inline
int
bicgstab( unsigned N, const MATRIX &A,
	  const double *b, double *x, double eps ) {
  return bicgstab(N,A,b,x,eps,false);
}

// ============================================================================


template< class MATRIX, class PC_MATRIX > inline
int
bicgstab( unsigned N, const MATRIX &A, const PC_MATRIX &C,
	  const double *b, double *x, double eps, bool detailed ) {
  double *rT  = new double[N];
  double *d   = new double[N];
  double *h   = new double[N];
  double *u   = new double[N];
  double *Sd  = new double[N];
  double *t   = new double[N];
  double *aux = new double[N];
  double *s   = h;
  double rTh, rTSd, rTr, alpha, beta, omega, st, tt;
  int its=0;
  double err=eps*eps*ddot(N,b,1,b,1);
  // f"ur's Abbruchkriterium (*)  -- r enth"alt immer das Residuum r=Ax-b
  double *r = new double[N];

  mult(A,x,r);
  daxpy(N,-1.,b,1,r,1);
  mult(C,r,d);
  dcopy(N,d,1,h,1);
  dcopy(N,h,1,rT,1);
  assert( ddot(N,rT,1,rT,1)>1e-40 );
  rTh=ddot(N,rT,1,h,1);
  rTr=ddot(N,r,1,r,1);
  while ( rTr>err ) {
    mult(A,d,aux);
    mult(C,aux,Sd);
    rTSd=ddot(N,rT,1,Sd,1);
    assert( fabs(rTSd)>1e-40 );
    alpha=rTh/rTSd;
    daxpy(N,-alpha,aux,1,r,1);
    dcopy(N,h,1,s,1);
    daxpy(N,-alpha,Sd,1,s,1);
    mult(A,s,aux);
    mult(C,aux,t);
    daxpy(N,1.,t,1,u,1);
    dscal(N,alpha,u,1);
    st=ddot(N,s,1,t,1);
    tt=ddot(N,t,1,t,1);
    if ( fabs(st)<1e-40 || fabs(tt)<1e-40 )
      omega = 0.;
    else
      omega = st/tt;
    daxpy(N,-omega,aux,1,r,1);
    daxpy(N,-alpha,d,1,x,1);
    daxpy(N,-omega,s,1,x,1);
    dcopy(N,s,1,h,1);
    daxpy(N,-omega,t,1,h,1);
    beta=(alpha/omega)/rTh; rTh=ddot(N,rT,1,h,1); beta*=rTh;
    dscal(N,beta,d,1);
    daxpy(N,1.,h,1,d,1);
    daxpy(N,-beta*omega,Sd,1,d,1);
    rTr=ddot(N,r,1,r,1);
    if ( detailed )
      cout<<"bicgstab "<<its<<"\t"<<sqrt(rTr)<<endl;
    ++its;
  }
  delete[] r;
  delete[] rT;
  delete[] d;
  delete[] h;
  delete[] u;
  delete[] Sd;
  delete[] t;
  delete[] aux;
  return its;
}

template< class MATRIX, class PC_MATRIX > inline
int
bicgstab( unsigned N, const MATRIX &A, const PC_MATRIX &C,
	  const double *b, double *x, double eps ) {
  return bicgstab(N,A,C,b,x,eps,false);
}


// ============================================================================

#endif // BICGSTAB_BLAS_H
