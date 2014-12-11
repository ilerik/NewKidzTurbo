#ifndef TURBO_GEOMFUNCTIONS
#define TURBO_GEOMFUNCTIONS

//Functions operating geometrical objects

#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "basetypes.h"
#include "cgnslib.h"
#include "optimization.h"
#include "mkl.h"
#include "mkl_lapack.h"

//Compute cgns element measure ( volume, square, length ) 
double ComputeElementMeasure() {
	return 0;
};

//Compute gradient using least squares
Vector ComputeGradientByPoints(Vector point, double value, const std::vector<Vector>& points, const std::vector<double>& values) {
		Vector grad;		

		//Build matrix (depends only on grid)
		/*alglib::real_1d_array rhs;
		alglib::real_1d_array x;
		alglib::ae_int_t info;
		alglib::densesolverlsreport rep;
		alglib::real_2d_array matrix;*/
		
		//Input
		int n = 3; //Number of dimensions (unknowns)
		int nPoints = points.size(); 
		int m = nPoints; //Number of equations
		int nrhs = 1; //Number of right hand side
		std::vector<double> a(n*m, 0);
		std::vector<double> b(nrhs*m, 0);
		std::vector<double> work(1 , 0);
		int lda = m; //lda = n;
		int ldb = m; //ldb = nrhs;
		std::vector<int> jpvt(n, 0);
		double rcond = 0.01;
		int lwork = -1;

		//Output
		int _info;
		int rank;				
				
		//matrix.setlength(points.size(), 3);
		for (int i = 0; i<points.size(); i++) {			
			Vector dr = point - points[i];
			/*matrix[i][0] = dr.x;
			matrix[i][1] = dr.y;
			matrix[i][2] = dr.z;*/

			a[i + 0*lda] = dr.x;
			a[i + 1*lda] = dr.y;
			a[i + 2*lda] = dr.z;
		};

		//rhs.setlength(points.size());		
		for (int i = 0; i<points.size(); i++) {						
			double dU = value - values[i];
			//rhs[i] = dU;

			b[0 * ldb + i] = dU;
		};

		//Workspace querry
		dgelsy(&m, &n, &nrhs, &a[0], &lda, &b[0], &ldb, &jpvt[0], &rcond, &rank, &work[0], &lwork, &_info);
		lwork = (int)work[0];
		work.resize(lwork, 0);

		//Solve problem
		dgelsy(&m, &n, &nrhs, &a[0], &lda, &b[0], &ldb, &jpvt[0], &rcond, &rank, &work[0], &lwork, &_info);
		if (_info != 0) throw Exception("Could not solve for gradient");
		grad = Vector(b[0], b[1], b[2]);
		//std::cout<<"grad = "<<grad.x<<" "<<grad.y<<" "<<grad.z<<"\n";

		/*x.setlength(3);
		alglib::rmatrixsolvels(matrix, points.size(), 3, rhs, 0.0, info, rep, x);
		if (info != 1) throw Exception("Could not solve for gradient");
		grad = Vector(x[0], x[1], x[2]);	*/
		//std::cout<<"grad_old = "<<grad.x<<" "<<grad.y<<" "<<grad.z<<"\n";

		//Output for check
	/*	std::cout<<"Matrix\n";
		for (int i = 0; i<points.size(); i++) {			
			Vector dr = point - points[i];
			std::cout<<dr.x<<" "<<dr.y<<" "<<dr.z<<"\n";			
		};
		std::cout<<"RHS\n";		
		for (int i = 0; i<points.size(); i++) {						
			double dU = value - values[i];
			std::cout<<dU<<"\n";	
		};*/

		return grad;
	};

//Compute gradient using least squares
Vector ComputeVelocityByPoints(int ndim, const std::vector<double>& velocities, const std::vector<Vector>& normals, const std::vector<double>& weights) {
	Vector velocity;			

	//Input
	int n = ndim; //Number of dimensions (unknowns)
	int nPoints = velocities.size(); 
	int m = nPoints; //Number of equations
	int nrhs = 1; //Number of right hand side
	std::vector<double> a(n*m, 0);
	std::vector<double> b(nrhs*m, 0);
	std::vector<double> work(1 , 0);
	int lda = m; //m; //lda = n;
	int ldb = m; //ldb = nrhs;
	std::vector<int> jpvt(n, 0);
	double rcond = 0.01;
	int lwork = -1;

	//Output
	int _info;
	int rank;				
				
	//Build matrix and right hand side
	for (int i = 0; i<velocities.size(); i++) {			
		Vector normal = normals[i];		
		double vn = velocities[i];
		b[0 * ldb + i]  = vn;
		if (ndim >= 1) a[i + 0*lda] = normal.x;
		if (ndim >= 2) a[i + 1*lda] = normal.y;
		if (ndim >= 3) a[i + 2*lda] = normal.z;
	};

	//Workspace querry
	dgelsy(&m, &n, &nrhs, &a[0], &lda, &b[0], &ldb, &jpvt[0], &rcond, &rank, &work[0], &lwork, &_info);
	lwork = (int)work[0];
	work.resize(lwork, 0);

	//Solve problem
	dgelsy(&m, &n, &nrhs, &a[0], &lda, &b[0], &ldb, &jpvt[0], &rcond, &rank, &work[0], &lwork, &_info);
	if (_info != 0) throw Exception("Could not solve for gradient");
	velocity = Vector(0, 0, 0);
	if (ndim >= 1) velocity.x = b[0];
	if (ndim >= 2) velocity.y = b[1];
	if (ndim >= 3) velocity.z = b[2];

	return velocity;
};

#endif