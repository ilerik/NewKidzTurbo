#ifndef TURBO_Utility_GeomFunctions
#define TURBO_Utility_GeomFunctions

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

//Determine angle between 2 vectors
double Angle(Vector r1, Vector r2) {
	double cosa = (r1 * r2) / (r1.mod() * r2.mod());
	return acos(cosa);
};

//Calc normal to polygon given by set of nodes
Vector CalcNormal(std::vector<Vector> points) {
	Vector norm;
	int np = points.size();	
	if (np == 2) {
		//2D case
		Vector r = points[0] - points[1];
		norm.x = r.y;
		norm.y = -r.x;
		norm /= norm.mod();
	} else {
		Vector center;
		//3D case
		for (int i = 0; i<np; i++) {
			center += points[i];
		};
		center /= np;
		for (int i = 0; i<np; i++) {
			Vector a = points[i] - center;
			Vector b = points[(i+1) % np] - center;
			norm += a & b;
		};
		norm /= norm.mod();
	}
	return norm;
};

//Calc rotation about axis on given angle
RotationMatrix CalcRotationAboutAxis(Vector v_axis, double angle) {
	RotationMatrix M;	
	v_axis /= v_axis.mod();	
	

	double angle_cos = cos(angle);
	double angle_sin = sin(angle);
	/*
	double angle_sin = (1 - angle * angle);	
	if (angle_sin > 0) { 
		angle_sin = sqrt(angle_sin);
	} else {
		angle_sin = 0;
	};
	if (angle_sin != angle_sin) throw new Exception("Bad angle_sin");
	*/
	M[0][0] = angle_cos + (1-angle_cos) * v_axis.x * v_axis.x;
	M[0][1] = (1-angle_cos) * v_axis.x * v_axis.y - angle_sin * v_axis.z;
	M[0][2] = (1-angle_cos) * v_axis.x * v_axis.z + angle_sin * v_axis.y;

	M[1][0] = (1-angle_cos) * v_axis.x * v_axis.y + angle_sin * v_axis.z;
	M[1][1] = angle_cos + (1-angle_cos) * v_axis.y * v_axis.y;	
	M[1][2] = (1-angle_cos) * v_axis.y * v_axis.z - angle_sin * v_axis.x;
	
	M[2][0] = (1-angle_cos) * v_axis.x * v_axis.z - angle_sin * v_axis.y;
	M[2][1] = (1-angle_cos) * v_axis.y * v_axis.z + angle_sin * v_axis.x;
	M[2][2] = angle_cos + (1-angle_cos) * v_axis.z * v_axis.z;
	return M;
};

//Calc rotation matrix to rotate one vector into another
RotationMatrix CalcRotation(const Vector& v1, const Vector& v2) {
	Vector v_axis = (v1 & v2);
	RotationMatrix M;
	if (v_axis.mod() < 1e-5) return M;
	double angle = (v1 * v2) / (v1.mod() * v2.mod());	//cos		
	if (angle != angle) {
		printf("%lg %lg\n", v1.mod(), v2.mod());
		printf("%lg %lg %lg\n", v2.x, v2.y, v2.z);
		throw Exception("Bad angle");
	};
	angle = acos(angle);
	M = CalcRotationAboutAxis(v_axis, angle);
	return M;
};

//Find intersection point of line & plane
Vector PlaneLineIntersec(const Vector P_norm, const Vector P_point, const Vector L1, const Vector L2)
{
	double k = P_norm * (P_point - L1);
	k /= P_norm * (L2 - L1);
	return L1 + k*(L2 - L1);
};

//Find Center of Points
Vector FindCenter(std::vector<Vector> points)
{
	Vector res = Vector(0, 0, 0);
	for(int i=0; i<points.size(); i++) res += points[i];
	res /= points.size();

	return res;
};

#endif