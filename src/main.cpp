#include "stdafx.h"

#include <stdio.h>
#include <string>
#include <sstream>
#include <iomanip>
#include "kernel.h"
#include "gridgeneration.h"
#include <cmath>
#include <random>
#include "MeshMovement.h"

#include "cuda.h"
#include <cuda_runtime.h>

//Test cases
#include "test_list.h"


template< typename T >
std::string int_to_hex( T i )
{
  std::stringstream stream;
  stream << "0x" 
         << std::setfill ('0') << std::setw(sizeof(T)*2) 
         << std::hex << i;
  return stream.str();
}


//Shear layer computation
class ImpactShockInitialConditions : public InitialConditions::InitialConditions
{
private:
	//State of target	
	double _u0;
	double _roL;
	double _roR;
	int _nmatLeft;
	int _nmatRight;
	//Speed of projectile	
	double _V;		
public:	
	ImpactShockInitialConditions(double V, int nmatLeft, double roL, int nmatRight, double roR) {
		_V = V;
		_roL = roL;
		_nmatLeft = nmatLeft;
		_roR = roR;
		_nmatRight = nmatRight;
	};

	virtual int getInitialGasModelIndex(const Cell& cell) {
		//Cell center
		double x = cell.CellCenter.x;
		double y = cell.CellCenter.y;
		double z = cell.CellCenter.z;
		
		if (x <= 0) {
			return _nmatLeft; //Pb
		} else {
			return _nmatRight; //Steel			
		};
	};

	virtual std::vector<double> getInitialValues(const Cell& cell) {
		int nmat = getInitialGasModelIndex(cell); //get material index

		std::vector<double> initValues;				
		//Other velocities
		double v = 0;
		double w = 0;

		//Internal energy
		double e = 0;

		//Left state		
		double roL = _roL;		
		double uL = 1.0*_V;
			
		//Right state		
		double roR = _roR;		
		double uR = -1.0*_V;

		//Cell center
		double x = cell.CellCenter.x;
		double y = cell.CellCenter.y;
		double z = cell.CellCenter.z;
				
		//Values
		double u = 0;
		double ro = 0;
		double roE = 0;				
		if (x <= 0) {
			ro = roL;
			u = uL;			
		} else {
			ro = roR;
			u = uR;						
		};
			
		//Convert to conservative variables
		//double P = 1e11;
		//e = P / (0.4 * ro);
		roE = ro*(e + (u*u + v*v + w*w) / 2.0);
		initValues.resize(_gasModels[nmat]->nConservativeVariables);
		initValues[0] = ro;
		initValues[1] = ro * u;
		initValues[2] = ro * v;
		initValues[3] = ro * w;
		initValues[4] = roE;// + (u*u + v*v + w*w) / 2.0); //TO DO check
		return initValues;
	};
};

void runImpactShockTest(int argc, char *argv[]) {
	Kernel _kernel;
	const double PI = 3.14159265359;
	
	_kernel.Initilize(&argc, &argv);
	double L = 15e-3; // 15 mm; 
	Grid _grid = GenGrid1D(_kernel.getParallelHelper(), 200, -L, L, false); //Change grid size here
	_kernel.BindGrid(&_grid);
	//_kernel.LoadGrid("C:\\Users\\Ilya\\Downloads\\cilindr5.11(1).cgns");

	//Fill in configuration
	_kernel.ReadConfiguration(""); //Change run parameters here
	_kernel.InitCalculation();

	//Initial conditions
	//double ro0 = 1;
	double roPb = 1000 * 1.0 / 0.88200003E-01; // SI	for Pb
	double roSteel = 1000 * 1.0 / 0.127; // SI	for stainless steel	
	double V = 250; //m/s
	///ImpactShockInitialConditions ic(V, 0, roPb, 0, roPb);
	//ImpactShockInitialConditions ic(V, 1, roSteel, 2, roPb);
	//ImpactShockInitialConditions ic(V, 1, roSteel, 1, roSteel);
	ImpactShockInitialConditions ic(V, 2, roPb, 2, roPb);
	_kernel.GenerateInitialConditions(&ic);	

	//Run test
	_kernel.RunCalculation();
	_kernel.FinalizeCalculation();

	//Output result
	_kernel.SaveGrid("result.cgns");
	_kernel.SaveSolution("result.cgns", "Solution");
	_kernel.Finalize();	
};

//void runImpactShockTest2D(int argc, char *argv[]) {
//	Kernel _kernel;
//	const double PI = 3.14159265359;
//	
//	_kernel.Initilize(&argc, &argv);
//	double L = 15e-3; // 15 mm; 
//	double xMin = -3*L;
//	double xMax = L;
//	double yMin = -L;
//	double yMax = L;
//	Grid _grid = GenGrid2D(_kernel.getParallelHelper(), 200, 100, xMin, xMax, yMin, yMax, 1.0, 1.0, false, false); //Change grid size here
//	//Modify grid for initial disturbances
//	double A = 1e-4;
//	std::default_random_engine generator;
//	std::uniform_real_distribution<double> distribution(-A,A);
//	std::vector<int> nodes;
//	std::vector<Vector> displacements;
//	for (Node& n : _grid.localNodes) {
//		//Distortion
//		if ((n.P.x == 0) && (std::abs(n.P.y) != L)) {
//			double delta = distribution(generator);
//			Vector dr(delta, 0, 0);
//			nodes.push_back(n.GlobalIndex);
//			displacements.push_back(dr);
//		};
//
//		//Unmovable borders
//		if ((n.P.x == xMin) || (n.P.x == xMax) || (n.P.y== yMin) || (n.P.y == yMax)) {
//			nodes.push_back(n.GlobalIndex);
//			displacements.push_back(Vector(0,0,0));
//		};
//	};
//
//	_kernel.BindGrid(&_grid);
//	_kernel.GenerateGridGeometry();
//	MeshMovement _moveHelper;
//	_moveHelper.IDWMove(_grid, nodes, displacements); 
//	//_kernel.LoadGrid("C:\\Users\\Ilya\\Downloads\\cilindr5.11(1).cgns");
//
//	//Fill in configuration
//	_kernel.ReadConfiguration(""); //Change run parameters here
//	_kernel.InitCalculation();
//
//	//Initial conditions
//	//double ro0 = 1;
//	double roPb = 1000 * 1.0 / 0.88200003E-01; // SI	for Pb
//	double roSteel = 1000 * 1.0 / 0.127; // SI	for stainless steel	
//	double V = 500; //m/s
//	//ImpactShockInitialConditions ic(V, 0, roSteel, 0, roPb);
//	ImpactShockInitialConditions ic(V, 1, roSteel, 2, roPb);
//	//ImpactShockInitialConditions ic(V, 1, roSteel, 1, roSteel);
//	//ImpactShockInitialConditions ic(V, 2, roPb, 2, roPb);
//	_kernel.GenerateInitialConditions(&ic);	
//
//	//Run test
//	_kernel.RunCalculation();
//	_kernel.FinalizeCalculation();
//
//	//Output result
//	_kernel.SaveGrid("result.cgns");
//	_kernel.SaveSolution("result.cgns", "Solution");
//	_kernel.Finalize();	
//};

//Main program ))
int main(int argc, char *argv[]) {	
 //   LomonosovFortovGasModel* gasModel = new LomonosovFortovGasModel(Logger());

	//GasModelConfiguration conf;		
	//conf.GasModelName = "LomonosovFortovGasModel";
	//conf.SetPropertyValue("MaterialIndex", 1);
	//conf.SetPropertyValue("SpecificHeatVolume", 130); //From http://www.diracdelta.co.uk/science/source/s/p/specific%20heat%20capacity/source.html#.VMr8MP6sXQI
	//conf.SetPropertyValue("MeltingTemperature", 600.622); //From http://www.diracdelta.co.uk/science/source/m/e/melting%20point/source.html#.VMr8x_6sXQI
	//double roMetal = 1000 * 1.0 / 0.88200003E-01;; //SI lead (Pb)
	//gasModel->loadConfiguration(conf);

	////double V = 8.820E-02;
	////double E = -3.240E-08;
	//double V = 6.132E-02;
	//double E = 6.855E-01;

	//GasModel::ConservativeVariables U;
	//U.ro = 1000.0 / V;//0.88200003E-01;
	//U.rou = 0;
	//U.rov = 0;
	//U.row = 0;
	//U.roE = U.ro * E * 1e6;
	//double C = 0;
	//double Gr= 0;
	//double P = 0; //gasModel->GetPressure(U);
	//double e = 0;
	//bool NF = false;
	//gasModel->EOSE5(1.0 / V, E, P, C, Gr, NF);
	//gasModel->EOSP5(1.0 / V, P, e, C, Gr, NF);
	//GasModel::MediumPhase phase = gasModel->GetPhase(U);

	//return 0;

	//TestCasesMetalsImpact::TestCaseMetalsImpact_1D_SteelVSSteel test(500);
	//TestCasesMetalsImpact::TestCaseMetalsImpact_1D_PbVSPb test(1000);

	//2:1 width ratio
	//TestCasesMetalsImpact::MetalsImpact1DTestCase test( 1000,
	//	80, //snapshots
	//	30e-3, //
	//	15e-3, 
	//	8e-6, // time = 8 mks
	//	TestCasesMetalsImpact::MetalType::StainlessSteel,
	//	0.0,
	//	TestCasesMetalsImpact::MetalType::Plumbum,
	//	-750.0
	//	);

	//4:1 width ratio
	//TestCasesMetalsImpact::MetalsImpact1DTestCase test( 1000,
	//	80, //snapshots
	//	60e-3, //
	//	15e-3, 
	//	8e-6, // time = 8 mks
	//	TestCasesMetalsImpact::MetalType::StainlessSteel,
	//	0.0,
	//	TestCasesMetalsImpact::MetalType::Plumbum,
	//	-750.0
	//	);

	//3:1 width ratio
	//TestCasesMetalsImpact::MetalsImpact1DTestCase test( 2000,
	//	80, //snapshots
	//	45e-3, //
	//	15e-3, 
	//	8e-6, // time = 8 mks
	//	TestCasesMetalsImpact::MetalType::StainlessSteel,
	//	0.0,
	//	TestCasesMetalsImpact::MetalType::Plumbum,
	//	-750.0
	//	);

	////3:1 width ratio
	////PbWidth = 2x15
	//double widthPb = 30e-3;
	//TestCasesMetalsImpact::MetalsImpact1DTestCase test( 2000,
	//	200, //snapshots
	//	3*widthPb, //
	//	widthPb, 
	//	20e-6, // time = 8 mks
	//	TestCasesMetalsImpact::MetalType::StainlessSteel,
	//	0.0,
	//	TestCasesMetalsImpact::MetalType::Plumbum,
	//	-750.0
	//	);

	//3:1 width ratio
	//PbWidth = 15, V = 1000
	//double widthPb = 15e-3;
	//TestCasesMetalsImpact::MetalsImpact1DTestCase test( 2000,
	//	100, //snapshots
	//	3*widthPb, //
	//	widthPb, 
	//	10e-6, // time = 8 mks
	//	TestCasesMetalsImpact::MetalType::StainlessSteel,
	//	0.0,
	//	TestCasesMetalsImpact::MetalType::Plumbum,
	//	-1000.0
	//	);

	//Test 1
	double widthSteel = 3e-3; //3 mm
	double widthPb = 2e-3; //2 mm
	//TestCasesMetalsImpact::MetalsImpact1DTestCase test( 500,
	//	300, //snapshots
	//	widthSteel, //
	//	widthPb, 
	//	3e-6, // time = 8 mks
	//	TestCasesMetalsImpact::MetalType::StainlessSteel,
	//	0.0,
	//	TestCasesMetalsImpact::MetalType::Plumbum,
	//	-470.0
	//	);

	//Test 2
	//widthSteel = 3e-3; //3 mm
	//widthPb = 9e-3; //2 mm
	//TestCasesMetalsImpact::MetalsImpact1DTestCase test( 1000,
	//	100, //snapshots
	//	widthSteel, //
	//	widthPb, 
	//	10e-6, // time = 8 mks
	//	TestCasesMetalsImpact::MetalType::StainlessSteel,
	//	0.0,
	//	TestCasesMetalsImpact::MetalType::Plumbum,
	//	-470.0
	//	);

	//Test 3
	//widthSteel = 30e-3; //3 mm
	//widthPb = 2e-3; //2 mm
	//TestCasesMetalsImpact::MetalsImpact1DTestCase test( 3000,
	//	100, //snapshots
	//	widthSteel, //
	//	widthPb, 
	//	10e-6, // time = 8 mks
	//	TestCasesMetalsImpact::MetalType::StainlessSteel,
	//	0.0,
	//	TestCasesMetalsImpact::MetalType::Plumbum,
	//	-470.0
	//	);

	//TestCases1D::TestCase1DALE3_RK4 test;
	//TestCases1D::TestCase1DALE4 test;

	//Test 1 2D
	double width = 1e-3;
	TestCasesMetalsImpact::MetalsImpact2DTestCase test( 10, 10,
		0, //snapshots
		width, //
		widthSteel, //
		widthPb, 
		3e-6, // time = 5 mks
		TestCasesMetalsImpact::MetalType::StainlessSteel,
		0.0,
		TestCasesMetalsImpact::MetalType::Plumbum,
		-470.0
		);
	test.RunTest(&argc, &argv);
	return 0;

};
