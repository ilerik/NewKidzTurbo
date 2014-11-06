#ifndef TURBO_BoundaryConditions_BoundaryCondition
#define TURBO_BoundaryConditions_BoundaryCondition

#include <vector>
#include <grid.h>
#include "GasModel.h"
#include "BoundaryConditionConfiguration.h"

namespace BoundaryConditions {

	//Boundary condition result type
	enum class BoundaryConditionResultType
	{
		Dirichlet, //Fixed value
		Neumann    //Fixed flux
	};

	//Base class for all boundary conditions
	class BoundaryCondition {
	public:	
		Grid* _grid;
		GasModel* _gasModel;	
	
		//Result types
		std::vector<BoundaryConditionResultType> bcResultTypes;

		virtual ~BoundaryCondition() {};

		//Interface functions
		virtual std::vector<double> getDummyValues(std::vector<double> values, const Cell& dummyCell) = 0;	

		virtual void setGrid(Grid& grid) {
			_grid = &grid;
		};

		virtual void setGasModel(GasModel& gasModel) {
			_gasModel = &gasModel;
		};

		//virtual void setRiemannSolver(RiemannSolver& rSolver) = 0;
		virtual void loadConfiguration(BoundaryConditionConfiguration& bcConfig) = 0; 
	};

};

#endif