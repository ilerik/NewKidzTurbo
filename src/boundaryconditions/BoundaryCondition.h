#ifndef TURBO_BoundaryConditions_BoundaryCondition
#define TURBO_BoundaryConditions_BoundaryCondition

#include <vector>
#include <grid.h>
#include "GasModel.h"
#include "BoundaryConditionConfiguration.h"

namespace BoundaryConditions {

	//Dummy values info
	struct DummyValuesInfo {
		std::vector<double> dummyValues;
		std::vector<double> dummyFluxes;
		int dummyMaterialIndex;
		int nVariables;
	};

	//Boundary condition result type
	enum class BoundaryConditionResultType
	{
		Dirichlet, //Fixed value
		Neumann    //Fixed gradient
	};	

	//Base class for all boundary conditions
	class BoundaryCondition {
	public:	
		int _nmat; //Gas model index
		Grid* _grid;
		std::vector<GasModel*> _gasModels;
	
		//Result types
		std::vector<BoundaryConditionResultType> bcResultTypes;

		//Boundary movement type
		BoundaryConditionMovementType movementType;

		virtual ~BoundaryCondition() {};

		//Interface functions
		/*virtual int getDummyGasModelIndex(const Cell& dummyCell) {			
			return nmat;
		};*/
		virtual std::vector<double> getDummyValues(int nmat, std::vector<double> values, const Cell& dummyCell) = 0;	

		virtual void setGrid(Grid& grid) {
			_grid = &grid;
		};

		virtual void setGasModel(std::vector<GasModel*>& gasModels) {
			_gasModels = gasModels;
		};

		virtual void setMaterialIndex(int nmat) {
			_nmat = nmat;
		};

		//Load configuration and check if it's applicable and correct
		virtual void loadConfiguration(BoundaryConditionConfiguration& bcConfig) = 0;
	};

};

#endif