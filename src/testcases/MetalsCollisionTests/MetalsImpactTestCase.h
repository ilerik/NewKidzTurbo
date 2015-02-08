#ifndef TURBO_TestCases_MetalsCollision_MetalsImpactTestCase
#define TURBO_TestCases_MetalsCollision_MetalsImpactTestCase

#include "TestCase.h"
#include "DSUClusteredSet.h"
#include "gengrid1D.h"

namespace TestCasesMetalsImpact {

//Enumeration of supported metal types
enum class MetalType {
	StainlessSteel,
	Plumbum
};

//Numerical test for ALE method from Luo et al 2004
//The well-known Sod shocktube problem is considered in this test case, whose solution contains 
//simultaneously a shock wave, a contact discontinuity, and an expansion fan.
class MetalsImpactTestCase: public TestCase {
protected:	

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

	BoundaryConditionConfiguration GetBoundaryConditionConfiguration(std::string materialName) {
		BoundaryConditionConfiguration conf;
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
				historyFile<<"VARIABLES = ";
				historyFile<<"\""<<"Time"<<"\" ";
				historyFile<<"\""<<"Iteration"<<"\" ";
				historyFile<<"\""<<"TimeStep"<<"\" ";
				historyFile<<"\""<<"MeltedZoneWidth"<<"\" ";
				historyFile<<"\""<<"TotalMeltedVolume"<<"\" ";
				historyFile<<"\""<<"xMin"<<"\" ";
				historyFile<<"\""<<"xMax"<<"\" ";
				historyFile<<"\""<<"xInterface"<<"\" ";
				historyFile<<"\""<<"xLeftBorder"<<"\" ";
				historyFile<<"\""<<"xRightBorder"<<"\" ";
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
			double xLeftBorder = 1000;
			double xRightBorder = -1000;
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
				historyFile<<xLeftBorder<<" ";
				historyFile<<xRightBorder<<" ";
				historyFile<<biggestClusterAvgTemperature<<" ";
				historyFile<<std::endl;
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
			if (x <= 0) {
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
			if (x <= 0) {
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