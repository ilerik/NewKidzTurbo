//  -*- C++ -*-

#ifndef GMRES_BLAS_H
#define GMRES_BLAS_H

// ============================================================================
//
//  GMRES nach Saad, Schultz
//     GMRES: a generalized minimal residual algorithm for solving nonsymmetric
//     linear systems
//     SIAM J Sci Stat Comput 7, 856-869 (1986)
//
//                                                 ----------------------------
//                                                 Christian Badura, Mai 1998
//
// ============================================================================


template< class Matrix >
inline int
gmres( int m, int N, const Matrix &A, const double *b, double *x, double eps );


template< class Matrix >
inline int
gmres( int m, int N, const Matrix &A, const double *b, double *x, double eps,
       bool detailed );


// ============================================================================


#include "mkl_blas.h"


template< class Matrix >
inline int
gmres( int m, int n, Matrix &A, const double *b, double *x, double eps,
       bool detailed ) {
  if ( n<=0 )
    return -1;

  typedef double *doubleP;
  double alpha = 0;
  const int one = 1;
  const char Transpose = 'T';
  const char NoTranspose = 'N';
  const char UpperTriangle = 'U';
  const char NotUnitTriangular = 'N';
  double *V  = new double[n*(m+1)];
  double *U  = new double[m*(m+1)/2];
  double *r  = new double[n];
  double *y  = new double[m+1];
  double *c  = new double[m];
  double *s  = new double[m];
  double **v = new doubleP[m+1];
  for ( int i=0; i<=m; ++i ) v[i]=V+i*n;
  int its=-1;
  {
    double beta, h, rd, dd, nrm2b;
    int j, io, uij, u0j;
    nrm2b=dnrm2(&n,b,&one);
    
    io=0;
    do  { // "aussere Iteration
      ++io;
      mult(A,x,r);
	  alpha = -1.0;
      daxpy(&n,&alpha,b,&one,r,&one);
      beta=dnrm2(&n,r,&one);
      dcopy(&n,r,&one,v[0],&one);
	  double a = 1./beta;
      dscal(&n,&a,v[0],&one);

      y[0]=beta;
      j=0;
      uij=0;
	do { // innere Iteration j=0,...,m-1
	u0j=uij;
	mult(A,v[j],v[j+1]);

	alpha = 1.0;
	beta = 0.0;
	int sy = j+1;
	dgemv(&Transpose,&n,&sy,&alpha,V,&n,v[j+1],&one,&beta,U+u0j,&one);

	alpha = -1.0;
	beta = 1.0;
	dgemv(&NoTranspose,&n,&sy,&alpha,V,&n,U+u0j,&one,&beta,v[j+1],&one);

	h=dnrm2(&n,v[j+1],&one);

	a = 1./h;
	dscal(&n,&a,v[j+1],&one);

	for ( int i=0; i<j; ++i ) { // rotiere neue Spalte
	  double tmp = c[i]*U[uij]-s[i]*U[uij+1];
	  U[uij+1]   = s[i]*U[uij]+c[i]*U[uij+1];
	  U[uij]     = tmp;
	  ++uij;
	}
	{ // berechne neue Rotation
	  rd     = U[uij];
	  dd     = sqrt(rd*rd+h*h);
	  c[j]   = rd/dd;
	  s[j]   = -h/dd;
	  U[uij] = dd;
	  ++uij;
	}
	{ // rotiere rechte Seite y (vorher: y[j+1]=0)
	  y[j+1] = s[j]*y[j];
	  y[j]   = c[j]*y[j];
	}
	++j;
	if ( detailed )
	  std::cout<<"gmres("<<m<<")\t"<<io<<"\t"<<j<<"\t"<<y[j]<<std::endl;
      } while ( j<m && fabs(y[j])>=eps*nrm2b );
      { // minimiere bzgl Y	
	int sx = j; 
	dtpsv(&UpperTriangle,&NoTranspose,&NotUnitTriangular,&sx,U,y,&one);
	// korrigiere X
	beta = 1.0;
	alpha = -1.0;
	dgemv(&NoTranspose,&n,&sx, &alpha,V, &n, y, &one, &beta, x, &one);
      }
    } while ( fabs(y[j])>=eps*nrm2b );

    // R"uckgabe: Zahl der inneren Iterationen
    its = m*(io-1)+j;
  }

  delete[] V;
  delete[] U;
  delete[] r;
  delete[] y;
  delete[] c;
  delete[] s;
  delete[] v;
  return its;
}


// ============================================================================


template< class Matrix >
inline int
gmres( int m, int n, const Matrix &A, const double *b, double *x, double eps ){
  return gmres(m,n,A,b,x,eps,false);
}

// ============================================================================


#endif // GMRES_BLAS_H
