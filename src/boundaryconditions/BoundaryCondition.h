#ifndef TURBO_GridLoading_CGNSReader
#define TURBO_BoundaryConditions_BoundaryCondition

#include <vector>

//Base class for all boundary conditions
class BoundaryCondition {
public:	
	virtual ~BoundaryCondition() {};

	//Interface functions
	virtual std::vector<double> getDummyCellValues() {} = 0;

private:

};

#endif