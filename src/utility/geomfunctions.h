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


//Compute cgns element measure ( volume, square, length ) 
double ComputeElementMeasure() {
	return 0;
};

//Compute gradient using least squares
Vector ComputeGradientByPoints(Vector point, double value, const std::vector<Vector>& points, const std::vector<double>& values) {
		Vector grad;		

		//Build matrix (depends only on grid)
		alglib::real_1d_array rhs;
		alglib::real_1d_array x;
		alglib::ae_int_t info;
		alglib::densesolverlsreport rep;
		alglib::real_2d_array matrix;
				
		matrix.setlength(points.size(), 3);
		for (int i = 0; i<points.size(); i++) {			
			Vector dr = point - points[i];
			matrix[i][0] = dr.x;
			matrix[i][1] = dr.y;
			matrix[i][2] = dr.z;
		};

		rhs.setlength(points.size());		
		for (int i = 0; i<points.size(); i++) {						
			double dU = value - values[i];
			rhs[i] = dU;
		};

		x.setlength(3);
		alglib::rmatrixsolvels(matrix, points.size(), 3, rhs, 0.0, info, rep, x);
		if (info != 1) throw Exception("Could not solve for gradient");
		grad = Vector(x[0], x[1], x[2]);	

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

#endif