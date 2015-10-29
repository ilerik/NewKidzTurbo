#ifndef TURBO_TestCases_MetalsCollision_MetalsImpactTestCase
#define TURBO_TestCases_MetalsCollision_MetalsImpactTestCase

#include "TestCase.h"
#include "DSUClusteredSet.h"
#include "gengrid1D.h"

namespace TestCasesMetalsImpact {

//Enumeration of supported metal types
enum class MetalType {
	StainlessSteel,
	Plumbum,
  Cuprum
};

//Numerical test for ALE method from Luo et al 2004
//The well-known Sod shocktube problem is considered in this test case, whose solution contains 
//simultaneously a shock wave, a contact discontinuity, and an expansion fan.
class MetalsImpactTestCase: public TestCase {
protected:	
	Grid _grid;					  //Grid object	
	Configuration _configuration; //Configuration object

	//Test parameters
	int _nSnapshots;
	double _widthLeft;
	double _widthRight;
	double _TimeMax;

	//Left state density, material index and velocity
	MetalType _metalLeft;
	double _uLeft;
	double _roLeft;
	int _nmatLeft;

	//Right state density, material index and velocity
	MetalType _metalRight;
	double _uRight;
	double _roRight;
	int _nmatRight;

	//Boundary condition (air)
	double _pBoundary;
	double _roBoundary;
	double _uBoundary;
	double _gammaBoundary;
	int _nmatBoundary;

	//Function for pertrubation position
	std::function<double(Vector)> getPertrubationPositionFunction;

public:
	//Default constructor
	MetalsImpactTestCase() {
	};

	//Parametrized constructor
	MetalsImpactTestCase(std::function<double(Vector)> pFunction, int nSnapshots, double widthLeft, double widthRight, double TimeMax, MetalType metalLeft, double uLeft, MetalType metalRight, double uRight) : 
		getPertrubationPositionFunction(pFunction),
		_nSnapshots(nSnapshots),
		_widthLeft(widthLeft),
		_widthRight(widthRight),
		_TimeMax(TimeMax),
		_metalLeft(metalLeft),
		_uLeft(uLeft),
		_metalRight(metalRight),
		_uRight(uRight)
	{
		_nmatBoundary = 2;
		_pBoundary = 1.1e5;
		_uBoundary = 0.0;
		_roBoundary = 1.225;
		_gammaBoundary = 1.4;
	};

	GasModelConfiguration GetGasModelConfiguration(MetalType metal, double &roMetal) {
		GasModelConfiguration conf;
		if (metal == MetalType::StainlessSteel) {
			conf.GasModelName = "LomonosovFortovGasModel";
			conf.SetPropertyValue("MaterialIndex", 0);
			conf.SetPropertyValue("SpecificHeatVolume", 510); //From http://www.diracdelta.co.uk/science/source/s/p/specific%20heat%20capacity/source.html#.VMr8MP6sXQI
			conf.SetPropertyValue("MeltingTemperature", 1800); //From http://www.diracdelta.co.uk/science/source/m/e/melting%20point/source.html#.VMr8x_6sXQI
			roMetal = 1000 * 1.0 / 0.127; //SI	for stainless steel;
			return conf;
		};
		if (metal == MetalType::Plumbum) {
			conf.GasModelName = "LomonosovFortovGasModel";
			conf.SetPropertyValue("MaterialIndex", 1);
			conf.SetPropertyValue("SpecificHeatVolume", 130); //From http://www.diracdelta.co.uk/science/source/s/p/specific%20heat%20capacity/source.html#.VMr8MP6sXQI
			conf.SetPropertyValue("MeltingTemperature", 600.622); //From http://www.diracdelta.co.uk/science/source/m/e/melting%20point/source.html#.VMr8x_6sXQI
			roMetal = 1000 * 1.0 / 0.88200003E-01;; //SI lead (Pb)
			return conf;
		};
    if (metal == MetalType::Cuprum) {
      conf.GasModelName = "LomonosovFortovGasModel";
      conf.SetPropertyValue("MaterialIndex", 2);
      conf.SetPropertyValue("SpecificHeatVolume", 130); //From http://www.diracdelta.co.uk/science/source/s/p/specific%20heat%20capacity/source.html#.VMr8MP6sXQI
      conf.SetPropertyValue("MeltingTemperature", 600.622); //From http://www.diracdelta.co.uk/science/source/m/e/melting%20point/source.html#.VMr8x_6sXQI
      roMetal = 1000 * 1.0 / 0.11200000E+00; //SI cuprum (Cu)
      return conf;
    };
	};

	GasModelConfiguration GetBoundaryGasModelConfiguration() {
		GasModelConfiguration conf;

		conf.GasModelName = "PerfectGasModel";
		conf.SetPropertyValue("IdealGasConstant", 8.3144621);
		conf.SetPropertyValue("SpecificHeatRatio", 1.4);
		conf.SetPropertyValue("SpecificHeatVolume", 1006.43 / 1.4);
		conf.SetPropertyValue("SpecificHeatPressure", 1006.43);
		
		return conf;		
	};

	enum class BoundaryConditionType {
		FreeSurface,
		FixedValues,
		Natural,
		Wall
	};

	BoundaryConditionConfiguration GetBoundaryConditionConfiguration(std::string materialName, BoundaryConditionType type) {
		BoundaryConditionConfiguration conf;
		if (type == BoundaryConditionType::FixedValues) {
			conf.BoundaryConditionType = BCType_t::BCInflowSupersonic;
			//conf.MovementType = BoundaryConditionMovementType::Fixed;
			conf.MovementType = BoundaryConditionMovementType::FreeSurface;
			conf.MaterialName = materialName;
			conf.SetPropertyValue("Density", _roBoundary);
			conf.SetPropertyValue("VelocityX", _uBoundary);
			conf.SetPropertyValue("VelocityY", 0);
			conf.SetPropertyValue("VelocityZ", 0);
			conf.SetPropertyValue("InternalEnergy", _pBoundary / ((_gammaBoundary - 1.0) * _roBoundary));
			return conf;
		};

		if (type == BoundaryConditionType::Wall) {
			conf.BoundaryConditionType = BCType_t::BCSymmetryPlane;
			conf.MovementType = BoundaryConditionMovementType::Fixed;
			conf.MaterialName = materialName;
			return conf;
		};

		if (type == BoundaryConditionType::FreeSurface) {			
			conf.BoundaryConditionType = BCType_t::BCGeneral;
			conf.MovementType = BoundaryConditionMovementType::FreeSurface;
			conf.MaterialName = materialName;
			return conf;
		};

		throw new Exception("Unspecified boundary condition type");
		return conf;
	};

	//Prepare computational grid
	virtual void PrepareGrid() = 0;

	//Prepare configuration
	virtual Configuration& PrepareConfiguration() = 0;	
		

	//History logger object
  class TestCaseHistoryLogger : public StepHistoryLogger {
  private:
    std::ofstream historyFile;    
  public:
    virtual void Init() {
      if (_parallelHelper->IsMaster()) {
        historyFile.open("history.dat", std::ofstream::out);
        historyFile << "VARIABLES = ";
        historyFile << "\"" << "Time" << "\" ";
        historyFile << "\"" << "Pressure, Cu (GPa)" << "\" ";
        historyFile << "\"" << "Pressure, Pb (GPa)" << "\" ";
        historyFile << "\"" << "Pressure, Interface (GPa)" << "\" ";
        historyFile << "\"" << "Velocity, Cu (m/s)" << "\" ";
        historyFile << "\"" << "Velocity, Pb (m/s)" << "\" ";
        historyFile << "\"" << "Velocity, Interface (m/s)" << "\" ";
        historyFile << "\"" << "Velocity, Free surface Cu (m/s)" << "\" ";        
        historyFile << "\"" << "Velocity, Free surface Pb (m/s)" << "\" ";
        historyFile << std::endl;
      };

      //Sync
      _parallelHelper->Barrier();
    };

    virtual void Finalize() {
      if (_parallelHelper->IsMaster()) {
        historyFile.close();
      };

      //Sync
      _parallelHelper->Barrier();
    };

    bool IsMelted(const int cellGlobalIndex) {
      int cellIndex = -1;
      int nmat = _kernel->GetCellGasModelIndex(cellGlobalIndex, cellIndex);
      GasModel::MediumPhase phase = _gasModels->at(nmat)->GetPhase(_kernel->GetCellValues(cellGlobalIndex, cellIndex));
      return phase == GasModel::MediumPhase::AboveMeltingPoint;
    };

    virtual void SaveHistory() {
      StepInfo* info = _kernel->getStepInfo();
      int iteration = info->Iteration;
      double time = info->Time;
      double timeStep = info->TimeStep;
   

      //Sensors   
      auto get_pressure = [&] (int cellIndex) {
        Cell* cell = _grid->localCells[cellIndex];
        int cellGlobalIndex = cell->GlobalIndex;
        double coordinate = cell->CellCenter.x;
        double volume = cell->CellVolume;
        int nmat = _kernel->GetCellGasModelIndex(cellGlobalIndex, cellIndex);
        double P = _gasModels->at(nmat)->GetPressure(_kernel->GetCellValues(cellGlobalIndex, cellIndex));
        double ro = (_kernel->GetCellValues(cellGlobalIndex, cellIndex))[0];
        double rou = (_kernel->GetCellValues(cellGlobalIndex, cellIndex))[1];
        double u = rou / ro;
        if (P < 0) P = 0; //Cutoff
        return P / 1e9; //GPa
      };

      auto get_velocity = [&] (int cellIndex) {
        Cell* cell = _grid->localCells[cellIndex];
        int cellGlobalIndex = cell->GlobalIndex;
        double coordinate = cell->CellCenter.x;
        double volume = cell->CellVolume;
        int nmat = _kernel->GetCellGasModelIndex(cellGlobalIndex, cellIndex);
        double P = _gasModels->at(nmat)->GetPressure(_kernel->GetCellValues(cellGlobalIndex, cellIndex));
        double ro = (_kernel->GetCellValues(cellGlobalIndex, cellIndex))[0];
        double rou = (_kernel->GetCellValues(cellGlobalIndex, cellIndex))[1];
        double u = rou / ro;
        return u; // m/s
      };

      const int NCells = this->_grid->nCells;
      const double widthCu = 3; // mm from surface
      const double xCu = 1.5; // mm from surface
      const double widthPb = 2; // mm from surface
      const double xPb = 1; // mm from surface
      const int iCu = std::floor(NCells * (widthPb + xCu) / (widthCu + widthPb));
      const int iPb = std::floor(NCells * (widthPb - xPb) / (widthCu + widthPb));
      const int iInterface = std::ceil(1 + NCells * (widthPb) / (widthCu + widthPb));

      double PCu{ get_pressure(iCu) };
      double PPb{ get_pressure(iPb) };
      double VCu{ get_velocity(iCu) };
      double VPb{ get_velocity(iPb) };
      double PInterface{ get_pressure(iInterface) };
      double VInterface{ get_velocity(iInterface) };
      double VFreeSurfacePb{ get_velocity(0) };      
      double VFreeSurfaceCu{ get_velocity(NCells - 1) };

      //Transformations      

      //Output
      if (_parallelHelper->IsMaster()) {
        historyFile << time * 1e6 << " ";
        historyFile << PCu << " ";
        historyFile << PPb << " ";
        historyFile << PInterface << " ";
        historyFile << VCu << " ";
        historyFile << VPb << " ";
        historyFile << VInterface << " ";
        historyFile << VFreeSurfaceCu << " ";
        historyFile << VFreeSurfacePb << " ";
        historyFile << std::endl;
      };

      //Sync
      _parallelHelper->Barrier();
    };

  };

	//Specify initial conditions
	class TestCaseInitialConditions : public InitialConditions::InitialConditions
	{
	private:		
	public:	
		//
		std::function<double(Vector)> _pFunction;

		//Left state density, material index and velocity
		double _uLeft;
		double _roLeft;
		int _nmatLeft;

		//Right state density, material index and velocity
		double _uRight;
		double _roRight;
		int _nmatRight;

		//Boundary states
		double _uBoundary;
		double _roBoundary;
		double _pBoundary;

		TestCaseInitialConditions(std::function<double(Vector)> pFunction, double uLeft, double roLeft, int nmatLeft, double uRight, double roRight, int nmatRight, double uBoundary, double roBoundary, double pBoundary) :
			_pFunction(pFunction),
			_uLeft(uLeft),
			_roLeft(roLeft),
			_nmatLeft(nmatLeft),
			_uRight(uRight),
			_roRight(roRight),
			_nmatRight(nmatRight),
			_uBoundary(uBoundary),
			_roBoundary(roBoundary),
			_pBoundary(pBoundary)
		{ };

		virtual int getInitialGasModelIndex(const Cell& cell) {
			//Cell center
			double x = cell.CellCenter.x;
			double y = cell.CellCenter.y;
			double z = cell.CellCenter.z;

			//Outside is air
			//if (x < -_widthLeft) return nmetAir;
			//if (x > LRight) return nmetAir;

			//Split in 2 halfs
			if ( _pFunction(cell.CellCenter) < 0) {
				return _nmatLeft; //First
			} else {
				return _nmatRight; //Second
			};		

			return 0;
		};

		virtual std::vector<double> getInitialValues(const Cell& cell) {
			int nmat = getInitialGasModelIndex(cell); //get material index

			std::vector<double> initValues;		
			double uL = _uLeft;
			double uR = _uRight;

			//Other velocities
			double v = 0;
			double w = 0;

			//Internal energy
			double e = 0;		

			//Cell center
			double x = cell.CellCenter.x;
			double y = cell.CellCenter.y;
			double z = cell.CellCenter.z;
				
			//Values
			double u = 0;
			double roE = 0;		
			double ro = 0;
			if ( _pFunction(cell.CellCenter) < 0) {
				u = uL;	
				ro = _roLeft;
			} else {
				u = uR;
				ro = _roRight;
			};

			//Outside is air
			/*double atm = 1.1e5;
			double gamma = 1.4;
			if ((x < - LLeft) || (x > LRight)) {
				ro = roAir;
				u = 0;
				e = atm / ((gamma - 1.0) * ro);
			};*/
			
			//Convert to conservative variables						
			roE = ro*(e + (u*u + v*v + w*w) / 2.0);
			initValues.resize(_gasModels[nmat]->nConservativeVariables);
			initValues[0] = ro;
			initValues[1] = ro * u;
			initValues[2] = ro * v;
			initValues[3] = ro * w;
			initValues[4] = roE;

			return initValues;
		};
	}; //Initial conditions

}; //TestCase

}; //namespace


#endif