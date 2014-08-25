#ifndef TURBO_InitialConditions_InitialConditions
#define TURBO_InitialConditions_InitialConditions

#include "grid.h"
#include "GasModel.h"

namespace InitialConditions {	

	//Base class for all custom initial conditions
	class InitialConditions {
		Grid* _grid;
		GasModel* _gasModel;	
	public:			
		virtual ~InitialConditions() {};

		//Interface functions
		virtual std::vector<double> getInitialValues(const Cell& cell) {
			std::vector<double> initValues;

			//For now 2D RTI (http://www.astro.virginia.edu/VITA/ATHENA/rt.html)
			//Grid sizes
			//double Lx = 0.02; // SI 2 cm
			//double Ly = 0.02; // SI 2 cm
			//double Lz = 0.02; // SI 2 cm
			double Lx = 0.5;
			double Ly = 1.5;

			//Acceleration
			//double a = -5e8; //SI
			double g = 1;			

			//Densities
			//double roDown = 7.9e3; //SI
			//double roUp = 11.34e3; //SI	
			double roDown = 1;
			double roUp = 2;
			double AtwoodNumber = (roUp - roDown)/(roDown+roUp);				

			//Cell center
			double x = cell.CellCenter.x;
			double y = cell.CellCenter.y;
			double z = cell.CellCenter.z;
			
			//Initial velocity pertrubation
			double A = 0.01;
			double u = 0;
			double v = A * (1+cos(2.*PI*x/Lx)) * (1+cos(2.*PI*y/Ly)) / 4.0; //1 mode
			double w = 0;

			//Initial pressure distribution 
			double Gamma = 1.4;			
			double P = 2.5;
			if (y <= 0) {
				P -= (roUp * g * y);
			} else {
				P -= (roDown * g * y);
			};
			//Specific stagnation energy
			double roE = P/(Gamma - 1.0);

			//Density distribution
			double ro;
			if (y <= 0) {
				ro = roDown;
			} else {
				ro = roUp;
			};
			
			//Convert to conservative variables
			initValues.resize(5);
			initValues[0] = ro;
			initValues[0] = ro * u;
			initValues[0] = ro * v;
			initValues[0] = ro * w;
			initValues[0] = roE;// + (u*u + v*v + w*w) / 2.0);

			return initValues;
		};	
		
		virtual void setGrid(Grid& grid) {
			_grid = &grid;
		};

		virtual void setGasModel(GasModel& gasModel) {
			_gasModel = &gasModel;
		};

		//virtual void setRiemannSolver(RiemannSolver& rSolver) = 0;
		virtual void loadConfiguration() {
		}; 
	};

};

#endif