#ifndef NewKidzTurbo_TestCases_RTInstabilityTests_RTInstabilityTestCase
#define NewKidzTurbo_TestCases_RTInstabilityTests_RTInstabilityTestCase

#include "TestCase.h"
#include "DSUClusteredSet.h"
#include "gridgeneration.h"

namespace RTInstabilityTests {

enum class GasModelType {
	LomonosovFortov,
	Barotropic,
	PerfectGas
};

enum class MaterialType {
	Heavy, //Hardcoded for now
	Light, //Hardcoded for now
	Helium,
	Xenon,
	StainlessSteel,
	Plumbum
};

enum class BoundaryConditionType {
	FreeSurface,
	Natural,
	Symmetry,
	Wall
};

struct GeometrySettings {
	//Mesh parameters
	int nCellX;
	int nCellY;

	//Geometry specific values
	double yMax;
	double yMin;
	double xMax;
	double xMin;

	//Disturbance parameters
	double yInterface;	//Position of undisturbed interface
	double A;			//Pertrubation amplitude
	double lambda;		//Wave length

	bool IsUserDefinedDisturbance;
	std::function<double(double)> pFunction;
};

//Materials settings
struct MaterialSettigs {
	GasModelType gasModelHeavy;
	double roHeavy; //
	GasModelType gasModelLight;
	double roLight; //
};

//Solution method settings
struct MethodSettings {
	std::string meshMotionType;	
	SpatialDiscretisationType spatialReconstruction;
};

//Test parameters
struct TestSettings {
	double Pinterface;	//Pressure at interface
	Vector gravity;		//Gravity field intensity
	BoundaryConditionType bcLeft;
	BoundaryConditionType bcRight;
	BoundaryConditionType bcTop;
	BoundaryConditionType bcBottom;

	//Calculation parameters
	int MaxIteration;
	int SaveSolutionSnapshotIterations;
	double MaxTime;
	double SaveSolutionSnapshotTime;

	//Other settings
	GeometrySettings geometrySettings;
	MaterialSettigs materialSettings;
	MethodSettings methodSettings;
};

//Numerical test for ALE method applie
//The well-known Sod shocktube problem is considered in this test case, whose solution contains 
//simultaneously a shock wave, a contact discontinuity, and an expansion fan.
class RTInstabilityTestCase: public TestCase,  public InitialConditions::InitialConditions {
protected:	
	std::shared_ptr<Grid> _gridPtr; //Grid object	
	Configuration _configuration; //Configuration object
	std::shared_ptr<ParallelManager> _MPIManager;  

	//Test settings
	TestSettings _settings;
	std::function<double(double)> _pFunction;
	std::function<double(Vector)> _signedDistanceFunction;

	Vector _interfaceCenter;
	Vector _interfaceNormal;
public:
	//Default constructor
	RTInstabilityTestCase() {
	};

	//Parametrized constructor
	RTInstabilityTestCase( ParallelManager& MPIManager, TestSettings settings )
	{
		_MPIManager = std::shared_ptr<ParallelManager>(&MPIManager);
		_settings = settings;

		//
		double xCenter = 0.5 * (_settings.geometrySettings.xMin + _settings.geometrySettings.xMax);
		_interfaceCenter = Vector(xCenter, _settings.geometrySettings.yInterface, 0);
		_interfaceNormal = Vector(0, 1.0, 0);

		//Initialize disturbance function
		_pFunction = [&](double x) {			
			double A = _settings.geometrySettings.A;
			double lambda = _settings.geometrySettings.lambda;
			double rI = A * std::cos(2 * PI * x / lambda);
			return rI;
		};

		//Signed distance function
		_signedDistanceFunction = [&](Vector r) {
			r = r - _interfaceCenter;
			double d = r * _interfaceNormal;
			Vector rSurface = r - d * _interfaceNormal;
			double x = rSurface.mod();
			double dSurface = _pFunction(x) + d;
			return dSurface;
		};

	};

	GasModelConfiguration GetGasModelConfiguration(GasModelType gasModel, MaterialType material) {
		GasModelConfiguration conf;
		if (gasModel == GasModelType::LomonosovFortov) {
			if (material == MaterialType::StainlessSteel) {
				conf.GasModelName = "LomonosovFortovGasModel";
				conf.SetPropertyValue("MaterialIndex", 0);
				conf.SetPropertyValue("SpecificHeatVolume", 510); //From http://www.diracdelta.co.uk/science/source/s/p/specific%20heat%20capacity/source.html#.VMr8MP6sXQI
				conf.SetPropertyValue("MeltingTemperature", 1800); //From http://www.diracdelta.co.uk/science/source/m/e/melting%20point/source.html#.VMr8x_6sXQI
				//roMetal = 1000 * 1.0 / 0.127; //SI	for stainless steel;
				return conf;
			};
			if (material == MaterialType::Plumbum) {
				conf.GasModelName = "LomonosovFortovGasModel";
				conf.SetPropertyValue("MaterialIndex", 1);
				conf.SetPropertyValue("SpecificHeatVolume", 130); //From http://www.diracdelta.co.uk/science/source/s/p/specific%20heat%20capacity/source.html#.VMr8MP6sXQI
				conf.SetPropertyValue("MeltingTemperature", 600.622); //From http://www.diracdelta.co.uk/science/source/m/e/melting%20point/source.html#.VMr8x_6sXQI
				//roMetal = 1000 * 1.0 / 0.88200003E-01;; //SI lead (Pb)
				return conf;
			};
		};

		if (gasModel == GasModelType::Barotropic) {
			double Pref = 2e14;	//P0 = 2*10^14 Pa	
			if (material == MaterialType::Light) {
				double roMetal = 4500; //steel reference density (ro2);
				conf.GasModelName = "BarotropicGasModel";
				conf.SetPropertyValue("YoungModulus", 2e13); // as for steel in http://en.wikipedia.org/wiki/Young%27s_modulus
				conf.SetPropertyValue("ReferencePressure", Pref); 
				conf.SetPropertyValue("ReferenceDensity", roMetal); 
				return conf;
			};
			if (material == MaterialType::Heavy) {
				double roMetal = 11400; //SI lead (Pb) ro1
				conf.GasModelName = "BarotropicGasModel";
				conf.SetPropertyValue("YoungModulus", 2e13); // as in http://en.wikipedia.org/wiki/Lead
				conf.SetPropertyValue("ReferencePressure", Pref); //
				conf.SetPropertyValue("ReferenceDensity", roMetal); //				
				return conf;
			};
		};

		if (gasModel == GasModelType::PerfectGas) {
			conf.GasModelName = "PerfectGasModel";
			conf.SetPropertyValue("IdealGasConstant", 8.3144621);
			double Pref = 1e10;	//P0 = 2*10^14 Pa
			double gammaHeavy = 1.4; // 
			double gammaLight = 1.4; // 
			double Cp = 1006.43;
			if ((material == MaterialType::Heavy) || (material == MaterialType::Light)) {
				conf.SetPropertyValue("SpecificHeatRatio", gammaHeavy);
				conf.SetPropertyValue("SpecificHeatVolume", Cp / gammaHeavy);
				conf.SetPropertyValue("SpecificHeatPressure", Cp);
				return conf;
			};

			if ((material == MaterialType::Heavy) || (material == MaterialType::Light)) {
				conf.SetPropertyValue("SpecificHeatRatio", gammaLight);
				conf.SetPropertyValue("SpecificHeatVolume", Cp / gammaLight);
				conf.SetPropertyValue("SpecificHeatPressure", Cp);
				return conf;
			};
		};

		//Error incompatible settings
		throw new Exception("Incompatible gas model and material settings");
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

	BoundaryConditionConfiguration GetBoundaryConditionConfiguration(std::string materialName, BoundaryConditionType type) {
		BoundaryConditionConfiguration conf;

		if (type == BoundaryConditionType::Wall) {
			conf.BoundaryConditionType = BCType_t::BCSymmetryPlane;
			conf.MovementType = BoundaryConditionMovementType::Fixed;
			conf.MaterialName = materialName;
			return conf;
		};

		if (type == BoundaryConditionType::Natural) {
			conf.BoundaryConditionType = BCType_t::BCGeneral;
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

	//Projects vector on interface and returns distance to it and interface reference frame coordinates of projection
	void GetDistanceAndInterfaceCoordinates(Vector& r, double& distance, std::vector<double>& coords) {
		Vector dr = (r - _interfaceCenter);
		Vector dn = (dr * _interfaceNormal) * _interfaceNormal;
		Vector dt = (dr - dn);
		distance = dr * _interfaceNormal;
		coords.clear();
		coords.push_back(dt.mod());
	};

	//Prepare computational grid
	void PrepareGrid() {
		_gridPtr = GenGrid2DInterfacePertrubation(_MPIManager, 
			_settings.geometrySettings.nCellX, _settings.geometrySettings.nCellY,
			_settings.geometrySettings.xMin, _settings.geometrySettings.xMax,
			_settings.geometrySettings.yMin, _settings.geometrySettings.yMax,
			1.0, 1.0,
			true, false,
			_interfaceCenter,
			_interfaceNormal,
			_signedDistanceFunction
			);

		_kernel->BindGrid(_gridPtr);
	};

	//Initial materials distribution
	virtual int getInitialGasModelIndex(const Cell& cell) {
		//Cell center
		double x = cell.CellCenter.x;
		double y = cell.CellCenter.y;
		double z = cell.CellCenter.z;

		//Split in 2 halfs
		Vector interfaceCenter = Vector(0.0, _settings.geometrySettings.yInterface, 0.0);
		Vector interfaceNormal = Vector(0.0, 1.0, 0.0);
		Vector dr = (cell.CellCenter - interfaceCenter);
		Vector dn = (dr * interfaceNormal) * interfaceNormal;
		double r = (dr - dn).mod();
		double signedDistance = _signedDistanceFunction(cell.CellCenter);
		if (signedDistance  < 0) {
			return 0; //Light
		} else {
			return 1; //Heavy
		};		

		return 0;
	};

	//Initial values distribution
	virtual std::vector<double> getInitialValues(const Cell& cell) {
		int nmat = getInitialGasModelIndex(cell); //get material index

		std::vector<double> initValues;		

		//Velocities
		double u = 0;
		double v = 0;
		double w = 0;

		//Internal energy
		double e = 0;		

		//Cell center
		double x = cell.CellCenter.x;
		double y = cell.CellCenter.y;
		double z = cell.CellCenter.z;
				
		//Values
		double roE = 0;		
		double ro = 0;
		double pressure = _settings.Pinterface;
		Vector dn = ((cell.CellCenter - _interfaceCenter) * _interfaceNormal) * _interfaceNormal;
		if (nmat == 0) {
			//Light
			ro = _settings.materialSettings.roLight;
			pressure += _settings.gravity * dn * ro;			
		} else {
			//Heavy
			ro = _settings.materialSettings.roHeavy;
			pressure += _settings.gravity * dn * ro;						
			v = -500;
		};
		e = pressure / (0.4 * ro);

		double Lx = _settings.geometrySettings.xMax - _settings.geometrySettings.xMin;
		double Ly = _settings.geometrySettings.yMax - _settings.geometrySettings.yMin;
		//v = 1e-4 * (std::cos(1200*x) * std::exp(-std::abs(y)));
			
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

	//Prepare configuration
	Configuration& PrepareConfiguration() {
		//Hardcode configuration for now
		_configuration.WorkingDirectory = "./results";
		_configuration.InputCGNSFile = "";
		_configuration.OutputCGNSFile = "result.cgns";

		//Rieman solver settings
		_configuration.RiemannSolverConfiguration.riemannSolverType = RiemannSolverConfiguration::RiemannSolverType::HLLC;

		//Availible gas models		
		//Left metal
		_configuration.AddGasModel("Light");
		_configuration.GasModelsConfiguration["Light"] = GetGasModelConfiguration(_settings.materialSettings.gasModelLight, MaterialType::Light);

		//Right metal
		_configuration.AddGasModel("Heavy");
		_configuration.GasModelsConfiguration["Heavy"] = GetGasModelConfiguration(_settings.materialSettings.gasModelHeavy, MaterialType::Heavy);

		//Boundary conditions		
		_configuration.BoundaryConditions["left"] = GetBoundaryConditionConfiguration("Light", RTInstabilityTests::BoundaryConditionType::Wall);
		_configuration.BoundaryConditions["right"] = GetBoundaryConditionConfiguration("Light", RTInstabilityTests::BoundaryConditionType::Wall);
		_configuration.BoundaryConditions["bottom"] = GetBoundaryConditionConfiguration("Light", _settings.bcBottom);
		_configuration.BoundaryConditions["top"] = GetBoundaryConditionConfiguration("Heavy", _settings.bcTop);
		
		//Solver settings					
		_configuration.SimulationType = TimeAccurate;
		_configuration.SpatialDiscretisation = _settings.methodSettings.spatialReconstruction;
		_configuration.CFL = 0.1;
		_configuration.RungeKuttaOrder = 4;		

		//ALE settings
		_configuration.ALEConfiguration.MeshMovementAlgorithm = MeshMovement::MeshMovementAlgorithm::IDW;
		_configuration.ALEConfiguration.ALEMotionType = _settings.methodSettings.meshMotionType;

		//Run settings
		_configuration.MaxIteration = _settings.MaxIteration;
		_configuration.MaxTime = _settings.MaxTime;
		_configuration.SaveSolutionSnapshotIterations = _settings.SaveSolutionSnapshotIterations;
		_configuration.SaveSolutionSnapshotTime = _settings.SaveSolutionSnapshotTime;

		//Gravity
		_configuration.g = _settings.gravity;

		//Bind configuration
		_kernel->VerboseOn();
		_kernel->BindConfiguration(_configuration);	

		return _configuration;
	};	

	//Get results of test run 
	virtual TestCaseResultInfo GetTestCaseResultInfo() override {
		TestCaseResultInfo info;
		throw new Exception("Not implemented");
		return info;
	};

	//Run kernel
	virtual void RunTestWithKernel(Kernel* kernel) override {
		//Set history logger
		//_kernel->setStepHistoryLogger(new TestCaseHistoryLogger());
		
		//Initialize calculation
		_kernel->InitCalculation();

		//Initial conditions
		_kernel->GenerateInitialConditions(this);
		_kernel->SaveGrid("init.cgns");
		_kernel->SaveSolution("init.cgns", "Solution");

		//Run computational cycle
		_kernel->RunCalculation();

		//Finilize
		_kernel->FinalizeCalculation();

		//Output result
		_kernel->SaveGrid("result.cgns");
		_kernel->SaveSolution("result.cgns", "Solution");
	};

	//Run test with program arguments
	virtual void RunTest(int* argc, char** argv[]) {
		_kernel = new Kernel();
		_kernel->Initilize(_MPIManager->comm());
		
		//Prepare grid and configuration
		this->PrepareGrid();
		this->PrepareConfiguration();

		//Call main function
		this->RunTestWithKernel(_kernel);

		_kernel->Finalize();
	}; 

	//History logger object
	/*class TestCaseHistoryLogger : public StepHistoryLogger {
	private:
		std::ofstream historyFile;
	public:
		virtual void Init() {
			if (_parallelHelper->IsMaster()) {
				historyFile.open("history.dat", std::ofstream::out);
				historyFile<<"VARIABLES = ";
				historyFile<<"\""<<"Time"<<"\" ";
				historyFile<<"\""<<"Iteration"<<"\" ";
				historyFile<<"\""<<"TimeStep"<<"\" ";
				historyFile<<"\""<<"MeltedZoneWidth"<<"\" ";
				historyFile<<"\""<<"TotalMeltedVolume"<<"\" ";
				historyFile<<"\""<<"xMin"<<"\" ";
				historyFile<<"\""<<"xMax"<<"\" ";
				historyFile<<"\""<<"xInterfaceMax"<<"\" ";
				historyFile<<"\""<<"xLeftBorder"<<"\" ";
				historyFile<<"\""<<"xRightBorder"<<"\" ";
				historyFile<<"\""<<"avgMeltedZoneTemperature"<<"\" ";
				historyFile<<"\""<<"xInterfaceMin"<<"\" ";
				historyFile<<"\""<<"dxInterface"<<"\" ";
				historyFile<<std::endl;
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
			StepInfo* info =_kernel->getStepInfo();
			int iteration = info->Iteration;
			double time = info->Time;
			double timeStep = info->TimeStep;
			double meltedZoneWidth = 0;
			double totalMeltedVolume = 0;
			double minPbCoordinate = 1000;
			double maxPbCoordinate = -1000;
			double minPbNotMeltedCoordinate = 1000;
			std::map<int, GasModel::MediumPhase> phases;
			std::map<int, int> nmats;
			std::vector<int> cellsIndexes;
			for (int cellIndex = 0; cellIndex < _gridPtr->nCellsLocal; cellIndex++) {
				Cell* cell = _gridPtr->localCells[cellIndex];
				int cellGlobalIndex = cell->GlobalIndex;
				cellsIndexes.push_back(cellGlobalIndex);
				double coordinate =  cell->CellCenter.x;
				double volume = cell->CellVolume;
				int nmat = _kernel->GetCellGasModelIndex(cellGlobalIndex, cellIndex);
				GasModel::MediumPhase phase = _gasModels->at(nmat)->GetPhase(_kernel->GetCellValues(cellGlobalIndex, cellIndex));
				phases[cellGlobalIndex] = phase;
				nmats[cellGlobalIndex] = nmat;

				//Total volume of fluid above melting point
				if (phase == GasModel::MediumPhase::AboveMeltingPoint) {				
					totalMeltedVolume += volume;
				};
			};

			//Initalize DSU structure
			DSUClusteredSet<int> meltedCells;			
			meltedCells.InitSet(cellsIndexes);

			//Iterate through all faces
			double xMin = 1000; bool isXMinSet = false;
			double xMax = 1000; bool isXMaxSet = false;
			bool isXInterfaceSet = false;
			double xInterfaceMin = 1000; 
			double xInterfaceMax = -1000;
			double xLeftBorder = 1000;
			double xRightBorder = -1000;
			for (Face& face : _gridPtr->localFaces) {
				bool isMeltedL = (phases[face.FaceCell_1] == GasModel::MediumPhase::AboveMeltingPoint);
				bool isMeltedR = (phases[face.FaceCell_2] == GasModel::MediumPhase::AboveMeltingPoint);
				bool isReversed = (face.FaceNormal.x < 0);
				int nmatL = nmats[face.FaceCell_1];
				int nmatR = nmats[face.FaceCell_2];
				double faceX = face.FaceCenter.x;
				if (_gridPtr->IsBoundaryFace(face)) continue;
				if (nmatL != nmatR) {
					if (xInterfaceMax < faceX) xInterfaceMax = faceX;
					if (xInterfaceMin > faceX) xInterfaceMin = faceX;
				};
				if (faceX > xRightBorder) xRightBorder = faceX;
				if (faceX < xLeftBorder) xLeftBorder = faceX;

				//Merge cells
				if ((isMeltedL && isMeltedR)) meltedCells.Union(face.FaceCell_1, face.FaceCell_2);
				if ((!isMeltedL && !isMeltedR)) meltedCells.Union(face.FaceCell_1, face.FaceCell_2);
			};

			//Process clusters
			std::function<bool(int)> ptr = std::bind(&TestCasesMetalsImpact::MetalsImpactTestCase::TestCaseHistoryLogger::IsMelted, this, std::placeholders::_1);
			std::vector< std::vector<int*> > clusters = meltedCells.GetClusters(ptr);
			int biggestCluster = -1;
			double biggestClusterLenght = 0;
			double biggestClusterAvgTemperature = 0;
			for (int i = 0; i < clusters.size(); i++) {
				double sumS = 0;
				double sumTS = 0;
				double xBegin = 1000;
				double xEnd = -1000;
				for (int j = 0; j < clusters[i].size(); j++) {
					int cellGlobalIndex = *clusters[i][j];
					int cellIndex = _gridPtr->cellsGlobalToLocal[cellGlobalIndex];
					Cell* cell = _gridPtr->localCells[cellIndex];
					double x =  cell->CellCenter.x;
					double S = cell->CellVolume;
					int nmat = nmats[cellGlobalIndex];
					double T = _gasModels->at(nmat)->GetTemperature(_kernel->GetCellValues(cellGlobalIndex, cellIndex));

					//Determine coordinates
					for (int faceIndex : cell->Faces) {
						Face& face = _gridPtr->localFaces[faceIndex];
						double faceX = face.FaceCenter.x;
						if (xBegin > faceX) xBegin = faceX;
						if (xEnd < faceX) xEnd = faceX;
					};

					//Average of temperature
					sumS += S;
					sumTS += T * S;
				};
				double avgT = sumTS / sumS;

				//Choose the biggest cluster
				if (sumS > biggestClusterLenght) {
					biggestClusterLenght = sumS;
					biggestClusterAvgTemperature = avgT;
					biggestCluster = i;
					xMax = xEnd;
					xMin = xBegin;
				};
			};

			if (biggestCluster == -1) {
				xMax = xMin = xInterfaceMax;
			};

			//Compute other quantities
			meltedZoneWidth = xMax - xMin;

			//Output
			if (_parallelHelper->IsMaster()) {
				historyFile<<time<<" ";
				historyFile<<iteration<<" ";
				historyFile<<timeStep<<" ";
				historyFile<<meltedZoneWidth<<" ";
				historyFile<<totalMeltedVolume<<" ";
				historyFile<<xMin<<" ";
				historyFile<<xMax<<" ";
				historyFile<<xInterfaceMax<<" ";
				historyFile<<xLeftBorder<<" ";
				historyFile<<xRightBorder<<" ";
				historyFile<<biggestClusterAvgTemperature<<" ";
				historyFile<<xInterfaceMin<<" ";
				historyFile<<xInterfaceMax - xInterfaceMin<<" ";
				historyFile<<std::endl;
			};

			//Sync
			_parallelHelper->Barrier();
		};
	}; */

}; //TestCase


/*//Specify initial conditions
class TestCaseInitialConditions : public InitialConditions::InitialConditions
{
private:		
	TestSettings& _settings;
	std::function<double(double)>& _pFunction;
public:	
	TestCaseInitialConditions( TestSettings& settings) :
		_settings(settings)
	{ 
		
	};

	virtual int getInitialGasModelIndex(const Cell& cell) {
		//Cell center
		double x = cell.CellCenter.x;
		double y = cell.CellCenter.y;
		double z = cell.CellCenter.z;

		//Split in 2 halfs
		Vector interfaceCenter;
		double r = (cell.CellCenter - interfaceCenter).mod();
		if ( _pFunction() < 0) {
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

		//Disturbance			
		u = 0;
		double Lx = 0.002;
		double Ly = 0.5e-3;
		double A = 0.1e-3;			
		double lambda = 2e-3;
		if ((x <= 0) && (std::abs(y) <= Ly / 2.0)) {
			double dL = std::abs(x);
			double unew = A * std::cos(2 * PI * y / lambda) * (1.0 - (dL / Lx));
			//if (unew > 0) unew = 0;
			u = unew;
		};
		//u = 0;
			
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
}; //Initial conditions */

}; //namespace

//Main

//RTInstabilityTests::TestSettings settings;

////Pertrubation
//settings.geometrySettings.lambda = 0.5;
//settings.geometrySettings.A = 0.0; //0.01; //0.1e-2;	
//settings.geometrySettings.IsUserDefinedDisturbance = false;

////Mesh
//int nWaves = 1;
//settings.geometrySettings.nCellX = 50;
//settings.geometrySettings.nCellY = 150;
//settings.geometrySettings.xMin = -0.5 * (nWaves * settings.geometrySettings.lambda);
//settings.geometrySettings.xMax = +0.5 * (nWaves * settings.geometrySettings.lambda);
//settings.geometrySettings.yMin = -0.75;
//settings.geometrySettings.yMax = +0.75;
//settings.geometrySettings.yInterface = 0.0;	

////Boundaries
//settings.bcTop = RTInstabilityTests::BoundaryConditionType::Wall;
//settings.bcBottom = RTInstabilityTests::BoundaryConditionType::Wall;

////Materials
//settings.materialSettings.gasModelHeavy = RTInstabilityTests::GasModelType::PerfectGas;
//settings.materialSettings.gasModelLight = RTInstabilityTests::GasModelType::PerfectGas;
//settings.materialSettings.roHeavy = 2.0;
//settings.materialSettings.roLight = 1.0;
//settings.Pinterface = 2.5;
////settings.gravity = Vector(0.0, 0.0, 0.0);
//settings.gravity = Vector(0.0, -0.25, 0.0);

////Method
//settings.methodSettings.meshMotionType = "ALEMaterialInterfaces";
////settings.methodSettings.meshMotionType = "Lagrangian";			
////settings.methodSettings.meshMotionType = "Eulerian";		
//settings.methodSettings.spatialReconstruction = SpatialDiscretisationType::WENO;

////Run parameters
//int nSnapshots = 100;
//settings.MaxTime = 10.0;
//settings.SaveSolutionSnapshotTime = settings.MaxTime / nSnapshots;
//settings.MaxIteration = 1000000;
//settings.SaveSolutionSnapshotIterations = 0;

////Create test object
//ParallelManager MPIManager(argc, argv);
//RTInstabilityTests::RTInstabilityTestCase test(MPIManager, settings);
//test.RunTest(&argc, &argv);


#endif