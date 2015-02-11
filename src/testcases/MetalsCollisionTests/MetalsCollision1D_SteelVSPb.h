#ifndef TURBO_TestCases_MetalsCollision_MetalsCollision1DSteelVSPb
#define TURBO_TestCases_MetalsCollision_MetalsCollision1DSteelVSPb

#include "TestCase.h"
#include "DSUClusteredSet.h"
#include "gengrid1D.h"

namespace TestCasesMetalsImpact {

//Numerical test for ALE method from Luo et al 2004
//The well-known Sod shocktube problem is considered in this test case, whose solution contains 
//simultaneously a shock wave, a contact discontinuity, and an expansion fan.
class TestCaseMetalsImpact_1D_SteelVSPb: public TestCase {
protected:
	Kernel* _kernel; //Computational kernel object
	Grid _grid;					  //Grid object	
	Configuration _configuration; //Configuration object
public:
	//Test parameters
	static const int nCells;
	static const double LLeft;
	static const double LRight;
	static const double TimeMax;

	//Left state density, material index for EOS5 and velocity
	static const double roSteel;
	static const int nmetSteel;
	static const double roPb;
	static const int nmetPb;
	static const double roAir;
	static const int nmetAir;


	//Impact relative speed
	double uImpact;

	//Parametrized constructor
	TestCaseMetalsImpact_1D_SteelVSPb(double _uImpact) {
		uImpact = _uImpact;
	};

	//Prepare computational grid
	void PrepareGrid() {		
		_grid = GenGrid1D(_kernel->getParallelHelper(), nCells, -LLeft, LRight, false);
		_kernel->BindGrid(&_grid);
	};

	//Prepare configuration object and set all parameters
	Configuration& PrepareConfiguration() {
		//Hardcode configuration for now
		_configuration.InputCGNSFile = "";
		_configuration.OutputCGNSFile = "result.cgns";

		//Rieman solver settings
		_configuration.RiemannSolverConfiguration.riemannSolverType = RiemannSolverConfiguration::RiemannSolverType::HLLC;

		//Availible gas models		

		//Stainless steel
		_configuration.AddGasModel("StainlessSteelLeft");
		_configuration.GasModelsConfiguration["StainlessSteelLeft"].GasModelName = "LomonosovFortovGasModel";
		_configuration.GasModelsConfiguration["StainlessSteelLeft"].SetPropertyValue("MaterialIndex", nmetSteel);
		_configuration.GasModelsConfiguration["StainlessSteelLeft"].SetPropertyValue("SpecificHeatVolume", 510); //From http://www.diracdelta.co.uk/science/source/s/p/specific%20heat%20capacity/source.html#.VMr8MP6sXQI
		_configuration.GasModelsConfiguration["StainlessSteelLeft"].SetPropertyValue("MeltingTemperature", 1800); //From http://www.diracdelta.co.uk/science/source/m/e/melting%20point/source.html#.VMr8x_6sXQI 

		_configuration.AddGasModel("PlumbumRight");
		_configuration.GasModelsConfiguration["PlumbumRight"].GasModelName = "LomonosovFortovGasModel";
		_configuration.GasModelsConfiguration["PlumbumRight"].SetPropertyValue("MaterialIndex", nmetPb);
		_configuration.GasModelsConfiguration["PlumbumRight"].SetPropertyValue("SpecificHeatVolume", 130); //From http://www.diracdelta.co.uk/science/source/s/p/specific%20heat%20capacity/source.html#.VMr8MP6sXQI
		_configuration.GasModelsConfiguration["PlumbumRight"].SetPropertyValue("MeltingTemperature", 600.622); //From http://www.diracdelta.co.uk/science/source/m/e/melting%20point/source.html#.VMr8x_6sXQI 

		//Air (ideal gas)
		double atm = 1.1e5;
		double gamma = 1.4;
		_configuration.AddGasModel("Air");
		_configuration.GasModelsConfiguration["Air"].GasModelName = "PerfectGasModel";
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("IdealGasConstant", 8.3144621);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatRatio", 1.4);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatVolume", 1006.43 / 1.4);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatPressure", 1006.43);

		//Boundary conditions
		//Left
		//_configuration.BoundaryConditions["left"].BoundaryConditionType = BCType_t::BCSymmetryPlane;
		_configuration.BoundaryConditions["left"].BoundaryConditionType = BCType_t::BCInflowSupersonic;
		//_configuration.BoundaryConditions["left"].MovementType = BoundaryConditionMovementType::Fixed;
		_configuration.BoundaryConditions["left"].MovementType = BoundaryConditionMovementType::FreeSurface;
		_configuration.BoundaryConditions["left"].MaterialName = "Air";
		_configuration.BoundaryConditions["left"].SetPropertyValue("Density", roAir);
		_configuration.BoundaryConditions["left"].SetPropertyValue("VelocityX", 0);
		_configuration.BoundaryConditions["left"].SetPropertyValue("VelocityY", 0);
		_configuration.BoundaryConditions["left"].SetPropertyValue("VelocityZ", 0);
		_configuration.BoundaryConditions["left"].SetPropertyValue("InternalEnergy", atm / ((gamma - 1.0) * roAir));

		//Right
		//_configuration.BoundaryConditions["right"].BoundaryConditionType = BCType_t::BCSymmetryPlane;
		_configuration.BoundaryConditions["right"].BoundaryConditionType = BCType_t::BCInflowSupersonic;
		//_configuration.BoundaryConditions["right"].MovementType = BoundaryConditionMovementType::Fixed;
		_configuration.BoundaryConditions["right"].MovementType = BoundaryConditionMovementType::FreeSurface;
		_configuration.BoundaryConditions["right"].MaterialName = "Air";
		_configuration.BoundaryConditions["right"].SetPropertyValue("Density", roAir);
		_configuration.BoundaryConditions["right"].SetPropertyValue("VelocityX", 0);
		_configuration.BoundaryConditions["right"].SetPropertyValue("VelocityY", 0);
		_configuration.BoundaryConditions["right"].SetPropertyValue("VelocityZ", 0);
		_configuration.BoundaryConditions["right"].SetPropertyValue("InternalEnergy", atm / ((gamma - 1.0) * roAir));
		
		//Solver settings					
		_configuration.SimulationType = TimeAccurate;
		_configuration.CFL = 0.5;
		_configuration.RungeKuttaOrder = 4;		

		//ALE settings
		//_configuration.ALEConfiguration.ALEMotionType = "Lagrangian";		
		_configuration.ALEConfiguration.ALEMotionType = "ALEMaterialInterfaces";		
		//_configuration.ALEConfiguration.ALEMotionType = "Eulerian";		

		//Run settings
		_configuration.MaxIteration = 1000000;
		_configuration.MaxTime = TimeMax;
		_configuration.SaveSolutionSnapshotIterations = 0;
		_configuration.SaveSolutionSnapshotTime = TimeMax / 80;			

		_kernel->BindConfiguration(_configuration);	

		return _configuration;
	};

	//Main interface function for running test case code
	void RunTestWithKernel(Kernel* kernel) {		
		_kernel = kernel; //Save reference to kernel
		
		//Get grid
		PrepareGrid();

		//Set parameters
		PrepareConfiguration();

		//Set history logger
		_kernel->setStepHistoryLogger(new TestCaseHistoryLogger());
		
		//Initialize calculation
		_kernel->InitCalculation();

		//Initial conditions
		_kernel->GenerateInitialConditions(new TestCaseInitialConditions(uImpact));	

		//Run computational cycle
		_kernel->RunCalculation();

		//Finilize
		_kernel->FinalizeCalculation();

		//Check results
		//Output result
		_kernel->SaveGrid("result.cgns");
		_kernel->SaveSolution("result.cgns", "Solution");
	};

	//Run test with program arguments
	void RunTest(int* argc, char** argv[]) {
		_kernel = new Kernel();
		_kernel->Initilize(argc, argv);

		//Call main function
		RunTestWithKernel(_kernel);

		_kernel->Finalize();
	};

	//History logger object
	class TestCaseHistoryLogger : public StepHistoryLogger {
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
				historyFile<<"\""<<"xInterface"<<"\" ";
				historyFile<<"\""<<"avgMeltedZoneTemperature"<<"\" ";
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
			for (int cellIndex = 0; cellIndex < _grid->nCellsLocal; cellIndex++) {
				Cell* cell = _grid->localCells[cellIndex];
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
			double xInterface = 1000; bool isXInterfaceSet = false;
			for (Face& face : _grid->localFaces) {
				bool isMeltedL = (phases[face.FaceCell_1] == GasModel::MediumPhase::AboveMeltingPoint);
				bool isMeltedR = (phases[face.FaceCell_2] == GasModel::MediumPhase::AboveMeltingPoint);
				bool isReversed = (face.FaceNormal.x < 0);
				int nmatL = nmats[face.FaceCell_1];
				int nmatR = nmats[face.FaceCell_2];
				double faceX = face.FaceCenter.x;
				if (_grid->IsBoundaryFace(face)) continue;
				if (nmatL != nmatR) {
					xInterface = faceX;
				};

				//Merge cells
				if ((isMeltedL && isMeltedR)) meltedCells.Union(face.FaceCell_1, face.FaceCell_2);
				if ((!isMeltedL && !isMeltedR)) meltedCells.Union(face.FaceCell_1, face.FaceCell_2);
			};

			//Process clusters
			std::function<bool(int)> ptr = std::bind(&TestCasesMetalsImpact::TestCaseMetalsImpact_1D_SteelVSPb::TestCaseHistoryLogger::IsMelted, this, std::placeholders::_1);
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
					int cellIndex = _grid->cellsGlobalToLocal[cellGlobalIndex];
					Cell* cell = _grid->localCells[cellIndex];
					double x =  cell->CellCenter.x;
					double S = cell->CellVolume;
					int nmat = nmats[cellGlobalIndex];
					double T = _gasModels->at(nmat)->GetTemperature(_kernel->GetCellValues(cellGlobalIndex, cellIndex));

					//Determine coordinates
					for (int faceIndex : cell->Faces) {
						Face& face = _grid->localFaces[faceIndex];
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
				xMax = xMin = xInterface;
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
				historyFile<<xInterface<<" ";
				historyFile<<biggestClusterAvgTemperature<<" ";
				historyFile<<std::endl;
			};

			//Sync
			_parallelHelper->Barrier();
		};
	};

	//Specify initial conditions
	//Sod's shock tube
	class TestCaseInitialConditions : public InitialConditions::InitialConditions
	{
	private:		
	public:	
		double uImpact;

		TestCaseInitialConditions(double _uImpact) {
			uImpact = _uImpact;
		};

		virtual int getInitialGasModelIndex(const Cell& cell) {
			//Cell center
			double x = cell.CellCenter.x;
			double y = cell.CellCenter.y;
			double z = cell.CellCenter.z;

			//Outside is air
			if (x < - LLeft) return nmetAir;
			if (x > LRight) return nmetAir;

			//Split in 2 halfs
			if (x <= 0) {
				return nmetSteel; //First
			} else {
				return nmetPb; //Second
			};		

			return 0;
		};

		virtual std::vector<double> getInitialValues(const Cell& cell) {
			int nmat = getInitialGasModelIndex(cell); //get material index

			std::vector<double> initValues;		
			double uL = 0;
			double uR = -uImpact;

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
			if (x <= 0) {
				u = uL;	
				ro = roSteel;
			} else {
				u = uR;
				ro = roPb;
			};

			//Outside is air
			double atm = 1.1e5;
			double gamma = 1.4;
			if ((x < - LLeft) || (x > LRight)) {
				ro = roAir;
				u = 0;
				e = atm / ((gamma - 1.0) * ro);
			};
			
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

//Test constant's
const int TestCaseMetalsImpact_1D_SteelVSPb::nCells = 1000;
const double TestCaseMetalsImpact_1D_SteelVSPb::LLeft = 15e-3; // 15 mm;
const double TestCaseMetalsImpact_1D_SteelVSPb::LRight = 15e-3; // 15 mm;
const double TestCaseMetalsImpact_1D_SteelVSPb::TimeMax = 8e-6; // 10 mks

//Density and material index for EOS5
const double TestCaseMetalsImpact_1D_SteelVSPb::roSteel = 1000 * 1.0 / 0.127; //SI	for stainless steel;
const int TestCaseMetalsImpact_1D_SteelVSPb::nmetSteel = 0;
const double TestCaseMetalsImpact_1D_SteelVSPb::roPb = 1000 * 1.0 / 0.88200003E-01;; //SI lead (Pb)
const int TestCaseMetalsImpact_1D_SteelVSPb::nmetPb = 1;
const double TestCaseMetalsImpact_1D_SteelVSPb::roAir = 1.225; //SI Air
const int TestCaseMetalsImpact_1D_SteelVSPb::nmetAir = 2;

}; //namespace


#endif