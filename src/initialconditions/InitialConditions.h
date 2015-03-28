#ifndef TURBO_InitialConditions_InitialConditions
#define TURBO_InitialConditions_InitialConditions

#include "grid.h"
#include "GasModel.h"
#include "configuration.h"

namespace InitialConditions {	

	//Base class for all custom initial conditions
class InitialConditions {		
public:			
	std::shared_ptr<Grid> _grid;
	std::vector<std::shared_ptr<GasModel>> _gasModels;

	virtual ~InitialConditions() {};		

	//Interface functions

	//Get material index number for cell
	virtual int getInitialGasModelIndex(const Cell& cell) = 0;

	//Generate initial values
	virtual std::vector<double> getInitialValues(const Cell& cell) = 0;
		
	virtual void setGrid(std::shared_ptr<Grid>& grid) {
		_grid = grid;
	};

	virtual void setGasModel(std::vector<std::shared_ptr<GasModel>>& gasModels) {
		_gasModels = gasModels;
	};

	//virtual void setRiemannSolver(RiemannSolver& rSolver) = 0;
	virtual void loadConfiguration() {
	}; 

};

};

#endif