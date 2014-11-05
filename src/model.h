#ifndef TURBO_MODEL
#define TURBO_MODEL

#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include "grid.h"
#include "datatypes.h"
#include "optimization.h"
#include "interpolation.h"
#include "geomfunctions.h"

//Step info
class StepInfo {
public:
	double TimeStep;
	std::vector<double> Residual;
	//maps from faces
	std::map<int, double> LyambdaC;
	std::map<int, double> LyambdaV;
	//maps from cells
	std::map<int, double> LyambdaSum;
};


//Medium properties
class MediumProperties {
public:
	double Gamma;
	double Cv;
	double Viscosity;
	double MolecularWeight;
	double ThermalConductivity;
};

//Wall information stored in each cell
struct CellWallInfo {
	int wallFaceIndex;
	double angle;
	double distance;
	double yPlus;
	double uPlus;
};

//Information stored in each wall face
struct FaceWallInfo {
	int wallFaceIndex;
	std::vector<int> layerCells;
	double boundaryLayerThickness;
	double shearVelocity;
	double shearStress;
};

//Identity that holds inside log part of boundary layer
void LogWallRelation(const alglib::real_1d_array &x, alglib::real_1d_array & f, void *ptr)
{
	std::vector<double>* params = (std::vector<double>*)ptr;
	double A = params->at(0);
	double B = params->at(1);
	//double shearVelocity = x[0];
	// target function
	/*if (shearVelocity < 1e-10) {
		f[0] = 1; 
	} else {
		f[0] = A + log(shearVelocity) - B / shearVelocity;
	};*/
	f[0] = A + x[0] - B/exp(x[0]);
}


//Class that represents model and its underlying properties
template<class RiemannSolver> 
class Model {
public:			
	//Boundary conditions section
	//Boundary condition types availible in model
	enum BCType{
		SubsonicInlet,
		SubsonicOutlet,
		NoSlip,
		Symmetry,
		Farfield
	};		

	class BoundaryCondition {		
		Dictionary parameters;
	public:
		Model& model;
		BCType type;
		
		void setParameter(std::string name, double value) {
			parameters.setValue(name, value);
		};	

		BoundaryCondition(Model& _model) : model(_model) {};
		virtual std::vector<double> ComputeConvectiveFlux(Face& f)
		{
			return model.ComputeConvectiveFlux(f);
		};

		virtual ConservativeVariables getDummyValues(ConservativeVariables UL, const Face& face) = 0;		
	};

	//Boundary conditions availible
	class NoSlipBoundaryCondition : public BoundaryCondition {
	public:
		NoSlipBoundaryCondition(Model& _model) : BoundaryCondition(_model) {};
		std::vector<double> ComputeConvectiveFlux(Face& f)	{
			ConservativeVariables UL = model.U[f.FaceCell_1];
			std::vector<double> res(5,0);
			double PL = (model.medium.Gamma - 1.0)*(UL.roE - 0.5*(UL.rou*UL.rou + UL.rov*UL.rov + UL.row*UL.row)/UL.ro);
			res[1] = PL*f.FaceNormal.x;
			res[2] = PL*f.FaceNormal.y;
			res[3] = PL*f.FaceNormal.z;

			return res;
		};

		ConservativeVariables getDummyValues(ConservativeVariables inV, const Face& face) {
			inV.rou *= -1;
			inV.rov *= -1;
			inV.row *= -1;
			return inV;
		};
	};

	class ConstantVelocityBoundaryCondition : public BoundaryCondition {
	public:
		Vector V;
		ConstantVelocityBoundaryCondition(Model& _model, Vector _V) : BoundaryCondition(_model), V(_V) {};

		ConservativeVariables getDummyValues(ConservativeVariables inV, const Face& face) {
			ConservativeVariables res = inV;
			Vector V_in(inV.rou, inV.rov, inV.row);
			V_in = (1.0/inV.ro)*V_in;

			Vector V_dummy = 2.0*V - V_in;
			res.rou = inV.ro*V_dummy.x;
			res.rov = inV.ro*V_dummy.y;
			res.row = inV.ro*V_dummy.z;
			res.roE = inV.roE - 0.5*inV.ro*V_in.mod()*V_in.mod() + 0.5*res.ro*V_dummy.mod()*V_dummy.mod();
			return res;
		};
	};

	class SymmetryBoundaryCondition : public BoundaryCondition {
	public:
		SymmetryBoundaryCondition(Model& _model) : BoundaryCondition(_model) {};

		ConservativeVariables getDummyValues(ConservativeVariables inV, const Face& face) {
			Vector roU(inV.rou, inV.rov, inV.row);			
			Vector newRoU = roU - 2*(roU * face.FaceNormal)*face.FaceNormal/(face.FaceNormal.mod()*face.FaceNormal.mod());
			inV.rou = newRoU.x;
			inV.rov = newRoU.y;
			inV.row = newRoU.z;
			return inV;
		};
	};

	class NaturalCondition : public BoundaryCondition {
	public:
		NaturalCondition(Model& _model) : BoundaryCondition(_model) {};

		ConservativeVariables getDummyValues(ConservativeVariables inV, const Face& face){
			return inV;
		};
	};

	//symmetry axis is the same that X axis
	class PeriodicAxisymmetricBoundaryCondition : public BoundaryCondition {
	public:
		PeriodicAxisymmetricBoundaryCondition(Model& _model) : BoundaryCondition(_model) {};
		double alpha;

		void setParams(double _alpha){
			alpha = _alpha;
		};

		ConservativeVariables getDummyValues(ConservativeVariables inV, const Face& face) {
			//first case
			ConservativeVariables res = inV;
			if(face.FaceCenter.z>0)
			{
				res.rov = inV.rov*cos(alpha) - inV.row*sin(alpha);
				res.row = inV.row*cos(alpha) + inV.rov*sin(alpha);
			}else
			{
				//second case
				res.rov = inV.rov*cos(alpha) + inV.row*sin(alpha);
				res.row = inV.row*cos(alpha) - inV.rov*sin(alpha);
			};
			return res;
		};
	};

	class SubsonicInletBoundaryCondition : public BoundaryCondition {
		double _pressure;
		double _temperature;

		////Velocity or velocity distribution
		Vector _velocity;
		//Velocity distribution
		Vector (*_velocityDistribution)( Vector, void * );	//Velocity distribution
		void *_velocityDistributionParams;					//Arbitrary parameters passed to main function

		MediumProperties medium;

	public:		
		SubsonicInletBoundaryCondition(Model& _model) : BoundaryCondition(_model) {			
			medium = model.medium;
			_velocityDistribution = NULL;
			_velocityDistributionParams = NULL;
		};

		void setParams(double pressure, double temperature, Vector velocity) {
			_pressure = pressure;
			_temperature = temperature;
			_velocity = velocity;
		};				

		void setVelocityDistribution( Vector (*vD)( Vector, void * ), void *params = NULL) {
			_velocityDistribution = vD;
			_velocityDistributionParams = params;
		};

		ConservativeVariables getDummyValues(ConservativeVariables inV, const Face& face) {			
			//Compute outgoing riemann invariant
			//Obtain speed of sound			
			Vector vd = model.GetVelocity(inV);
			double cd = model.GetSoundSpeed(inV);
			double Rminus = vd * face.FaceNormal / face.FaceNormal.mod() -  2*cd/(medium.Gamma - 1.0);
			//Now solve equation Rminus = const for velocity
			//Use Blazek recomendations	
			/*double cost = -1.0;
			if (vd.mod() > std::numeric_limits<double>::epsilon()) cost = - vd * face.FaceNormal / (face.FaceNormal.mod() * vd.mod());
			double cd2 = cd*cd;
			double c02 = cd2 + (medium.Gamma - 1.0) * vd * vd / 2.0;
			double A = (medium.Gamma - 1.0) * cost * cost + 2;
			double cb = A * c02 / (Rminus * Rminus * (medium.Gamma - 1.0)) - (medium.Gamma - 1.0) / 2.0;
			cb = sqrt(cb);
			cb = 1 + cost * cb;
			cb *= (-Rminus *  (medium.Gamma - 1.0)) / A;
			double cb2 = cb * cb;
			double T0 = _temperature;	
			double k = cb2/c02;
			if (k > 1.0) k = 1.0;
			double Tb = T0 * k;
			double P0 = _pressure;
			double Pb = P0 * pow(Tb/T0, medium.Gamma / (medium.Gamma - 1.0) );
			double vbmod = sqrt(2* medium.Gamma * medium.Cv * (T0 - Tb));
			Vector vb = _velocity * vbmod / _velocity.mod();*/		
			Vector vb = _velocity;
			if (_velocityDistribution != NULL) {
				vb = _velocityDistribution(face.FaceCenter, _velocityDistributionParams);
			};
			double P0 = _pressure;	//TO DO лишнии переменные
			double Pb = P0;	
			double Tb = _temperature;
			return model.PrimitiveToConservativeVariables(vb, Pb, Tb, medium);
		};
	};

	class SupersonicInletBoundaryCondition : public BoundaryCondition {
		double _pressure;
		double _temperature;
		Vector _velocity;		
		MediumProperties medium;
	public:		
		SupersonicInletBoundaryCondition(Model& _model) : BoundaryCondition(_model) {			
			medium = model.medium;
		};

		void setParams(double pressure, double temperature, Vector velocity) {
			_pressure = pressure;
			_temperature = temperature;
			_velocity = velocity;
		};			

		void setVelocityDistribution( Vector (*vD)( Vector, void * )) {
			_velocityDistribution = vD;
		};

		ConservativeVariables getDummyValues(ConservativeVariables inV, const Face& face) {												
			double Pb = _pressure;
			double Tb = _temperature;			
			Vector vb = _velocity;				
			return model.PrimitiveToConservativeVariables(vb, Pb, Tb, medium);
		};
	};

	class InletBoundaryCondition : public BoundaryCondition {
		double _pressure;
		double _temperature;
		Vector _velocity;		
		MediumProperties medium;

		//Boundary condition to choose
		SupersonicInletBoundaryCondition supersonicBC;
		SubsonicInletBoundaryCondition subsonicBC;
	public:		
		InletBoundaryCondition(Model& _model) : BoundaryCondition(_model), supersonicBC(_model), subsonicBC(_model) {			
			medium = model.medium;					
		};

		void setParams(double pressure, double temperature, Vector velocity) {
			_pressure = pressure;
			_temperature = temperature;
			_velocity = velocity;
			supersonicBC.setParams(pressure, temperature, velocity);
			subsonicBC.setParams(pressure, temperature, velocity);
		};						

		ConservativeVariables getDummyValues(ConservativeVariables inV, const Face& face) {												
			double cd = model.GetSoundSpeed(inV);
			Vector vd = model.GetVelocity(inV);
			if (cd < -vd * face.FaceNormal / face.FaceNormal.mod()) {
				return supersonicBC.getDummyValues(inV, face);
			} else {
				return subsonicBC.getDummyValues(inV, face);
			};
		};
	};

	class SubsonicOutletBoundaryCondition : public BoundaryCondition {
		double _pressure;
		MediumProperties medium;
	public:
		SubsonicOutletBoundaryCondition(Model& _model) : BoundaryCondition(_model) {			
			medium = model.medium;
		};

		void setParams(double pressure) {
			_pressure = pressure;
		};		

		ConservativeVariables getDummyValues(ConservativeVariables inV, const Face& face) {												
			//Now solve equation Rminus = const for velocity
			//Use Blazek recomendations						
			double Tb = model.GetTemperature(inV);			
			double Pb = _pressure;			
			Vector vb = model.GetVelocity(inV);			
			return model.PrimitiveToConservativeVariables(vb, Pb, Tb, medium);
		};
	};

//protected:
	//Internal variables
	Grid _grid;	//Grid
	std::map<int, BoundaryCondition*> _boundaryConditions; //Boundary conditions	

	//Store conservative variables for each cell
	DistributedEntityManager<ConservativeVariables> U;	

	//Riemann solver object
	RiemannSolver rSolver;	

	alglib::kdtree kdtWall; 	//KDTree structure and related operations for wall faces
	alglib::kdtree kdtCells; 	//KDTree structure and related operations for all cells

	std::map<int, CellWallInfo> _wallInfo;	//cellIndex -> wall information
	std::map<int, FaceWallInfo> _wallFaces; //wall faceIndex -> wall information
	std::map<int, double> _boundaryLayerHeight; //
	std::set<int> _wallBoundaryMarkers;
	inline bool IsWallFace(Face& face) {
		if ((face.isExternal) && ( _wallBoundaryMarkers.find(face.BCMarker) !=  _wallBoundaryMarkers.end() )) return true;
		return false;
	};

	//Computational settings
	double CFL;	//CFL condition value (0.35 for example)
	bool IsGradientRequired;
	bool IsViscousFluxesRequired;
	bool IsAxySymmetric;
	bool IsLimiterRequired;

	//Computed variables for each face
	std::map<int, std::vector<double>> fluxes;
	std::map<int, std::vector<double>> vfluxes; //TO DO Remove	

	//Gradient data storage
	std::map<int, std::set<int>> cellNeigbours;
	std::map<int, std::vector<Vector>> gradsFaces;		
	std::map<int, std::vector<Vector>> gradsCells;	
	std::map<int, Vector> gradCellsU;	
	std::map<int, Vector> gradFacesU;	
	std::map<int, Vector> gradCellsV;	
	std::map<int, Vector> gradFacesV;	
	std::map<int, Vector> gradCellsW;	
	std::map<int, Vector> gradFacesW;	
	std::map<int, Vector> gradCellsT;	
	std::map<int, Vector> gradFacesT;
	//Conservative variables gradients
	std::map<int, Vector> gradCellsRo;
	std::map<int, Vector> gradCellsRoU;
	std::map<int, Vector> gradCellsRoV;
	std::map<int, Vector> gradCellsRoW;
	std::map<int, Vector> gradCellsRoE;

	//limiters map
	std::map<int, std::vector<double> > limiters;		//from cell Index to limiters values
	
	//Residual information
	std::map<int, double> rouResidual;
	std::map<int, double> rovResidual;
	std::map<int, double> rowResidual;
	std::map<int, double> roResidual;
	std::map<int, double> roEResidual;

public:
	//Medium properties
	MediumProperties medium;
	double SoartingAngle;

	//external forces
	Vector g;	//uniform acceleration

	//Properties
	StepInfo stepInfo;
	double totalTime;
	int SchemeOrder;

	//Accessors and settings
	void SetGamma(double gamma) {
		medium.Gamma = gamma;
		rSolver.SetGamma(gamma);
	};

	void SetCv(double cv) {
		medium.Cv = cv;
	};

	void SetViscosity(double v) {
		medium.Viscosity = v;
	};

	void SetMolecularWeight(double m) {
		medium.MolecularWeight = m;
	};

	void SetThermalConductivity(double tk) {
		medium.ThermalConductivity = tk;
	};

	void SetCFLNumber(double cfl) {
		CFL = cfl;
	};

	void SetHartenEps(double eps)
	{
		rSolver.SetHartenEps(eps);
	};

	void SetUniformAcceleration(Vector _g) {
		g = _g;
	};

	//TO DO
	void EnableViscous() {
		IsGradientRequired = true;
		IsViscousFluxesRequired = true;
	};

	void DisableViscous() {
		IsGradientRequired = false;
		IsViscousFluxesRequired = false;
	};

	void EnableLimiter() {
		IsLimiterRequired = true;
	};

	void DisableLimiter() {
		IsLimiterRequired = false;
		std::vector<double> unity;						//fill vector of size 5 by unities
		for(int i=0; i<5; i++) unity.push_back(1);

		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		for(int i=0; i<cells.size(); i++)
		{
			limiters[cells[i]->GlobalIndex] = unity;
		};
		return;
	};

	void SetAngle(double angle)
	{
		SoartingAngle = angle;
	};

	void SetSchemeOrder(int order){
		SchemeOrder = order;
	};

	void SetAxiSymmetric()
	{
		IsAxySymmetric = true;
	};


	//Constructor
	Model() {	
		totalTime = 0; //Init simulation time
		IsAxySymmetric = false;	//default value
		IsLimiterRequired = false;
		g = Vector(0, 0, 0);
	};

	void BindGrid(Grid& grid) {
		//Grid
		_grid = grid;

		//Allocate memory for data over each cell		
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		for (int i = 0; i<cells.size(); i++) {					
			U.add(ConservativeVariables(cells[i]->GlobalIndex));
		};
	};

	void SetBoundaryCondition(std::string BCName, BoundaryCondition& bc) {
		if (_grid.patchesNames.find(BCName) == _grid.patchesNames.end()) throw Exception("No such boundary availible");
		int bcMarker = _grid.patchesNames[BCName];
		_boundaryConditions[bcMarker] = &bc; 
	};
	
	void SetInitialConditions(ConservativeVariables initValue) {
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		for (int i = 0; i<cells.size(); i++) {
			U[cells[i]->GlobalIndex] = initValue;
			U[cells[i]->GlobalIndex].GlobalIndex = cells[i]->GlobalIndex;
		};
		return;
	};
	void SetInitialConditions(ConservativeVariables(*funcInitValue)(Vector, void *), void *params = NULL) {
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		for (int i = 0; i<cells.size(); i++) {
			U[cells[i]->GlobalIndex] = funcInitValue(cells[i]->CellCenter, params);
			U[cells[i]->GlobalIndex].GlobalIndex = cells[i]->GlobalIndex;
		};
		return;
	};

	void SetInitialConditionsLinear(ConservativeVariables initValue) {
		const double dUdy = 5;
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();		
		for (int i = 0; i<cells.size(); i++) {
			ConservativeVariables cellValue;
			cellValue = initValue;
			cellValue.rou = initValue.ro * cells[i]->CellCenter.y * dUdy;
			U[cells[i]->GlobalIndex] = cellValue;
			U[cells[i]->GlobalIndex].GlobalIndex = cells[i]->GlobalIndex;
		};
		return;
	};

	void SetWallBoundary(std::string BCName, bool isWall) {
		//Obtain boundary bc marker
		if (_grid.patchesNames.find(BCName) == _grid.patchesNames.end()) throw Exception("No such boundary availible");
		int bcMarker = _grid.patchesNames[BCName];

		if (isWall) {
			_wallBoundaryMarkers.insert(bcMarker);
		
			//Move through all faces and add wall face if needed
			for (std::set<int>::iterator it = _grid.patches[bcMarker].faces_idx.begin(); it != _grid.patches[bcMarker].faces_idx.end(); it++) {
				_wallFaces[*it] = FaceWallInfo();
				_wallFaces[*it].wallFaceIndex = *it;
			};
		} else {
			if (_wallBoundaryMarkers.find(bcMarker) != _wallBoundaryMarkers.end()) _wallBoundaryMarkers.erase(bcMarker);

			//Move through all faces and remove wall face if needed
			for (std::set<int>::iterator it = _grid.patches[bcMarker].faces_idx.begin(); it != _grid.patches[bcMarker].faces_idx.end(); it++) {
				if (_wallFaces.find(*it) != _wallFaces.end()) _wallFaces.erase(_wallFaces.find(*it));
			};
		};
	};

	void Init() {	
		//ComputeFunctionGradient<ConservativeVariables>(gradCellsU, U, &Model<RiemannSolver>::GetVelocityX);
		std::cout<<GetPressure(U[1])<<"\n";
		std::cout<<GetVelocityX(U[1])<<"\n";
		std::cout<<GetTemperature(U[1])<<"\n";
	};

	//empty function for general case
	virtual void AditionalStepComputations()
	{
	};

	//Limiters Part
	std::vector<double> ComputeLimiter(Cell& c) {

	ConservativeVariables U_c = U[c.GlobalIndex];		//Values of c cell
	ConservativeVariables U_max = U_c;					//Max and Min Values from c cell and its neighbours
	ConservativeVariables U_min = U_c;			

	for(int i = 0; i<c.Faces.size(); i++)			//find U_min and U_max
	{
		ConservativeVariables U_neigh;						//Values in i-neighbour cell
		int neigh_cell;									//index of neighbour cell
		Face cur_face = _grid.faces[c.Faces[i]];
		if(cur_face.isExternal!=0) continue;
		if(cur_face.FaceCell_1 == c.GlobalIndex)		//find i-neighbour cell
		{
			neigh_cell = cur_face.FaceCell_2;
		}else neigh_cell = cur_face.FaceCell_1;

		U_neigh = U[neigh_cell];				//Values of i-neighbour cell
		for (int j = 0; j < 5; j++)
		{
				if(U_max[j] < U_neigh[j]) U_max[j] = U_neigh[j];	//find U_max,U_min
				if(U_min[j] > U_neigh[j]) U_min[j] = U_neigh[j];
		};
	};

	std::vector<double> res(5,1);			//Limiter values for all conservetive variables in Cell c					
	std::vector<Vector> grad;				//Gradient in Cell c
	grad.push_back(gradCellsRo[c.GlobalIndex]);
	grad.push_back(gradCellsRoU[c.GlobalIndex]);
	grad.push_back(gradCellsRoV[c.GlobalIndex]);
	grad.push_back(gradCellsRoW[c.GlobalIndex]);
	grad.push_back(gradCellsRoE[c.GlobalIndex]);

	for(int i=0; i<c.Faces.size(); i++)			//calc limiter - res
	{
		std::vector<double> res_new(5);			//to find minimal values
		Vector FaceCenter = _grid.faces[c.Faces[i]].FaceCenter;
		std::vector<double> delta2(5);							//special function					
		for(int j=0; j<5; j++)
		{
			delta2[j] = grad[j]*(FaceCenter - c.CellCenter);
			if(delta2[j] > 0) res_new[j] = (U_max[j] - U_c[j])/delta2[j];
			if(delta2[j] < 0) res_new[j] = (U_min[j] - U_c[j])/delta2[j];
			if(delta2[j] ==0) res_new[j] = 1;
			if(res_new[j] > 1)	res_new[j] = 1;
			if(res[j] > res_new[j]) res[j] = res_new[j];	//find min res_new
		};
	};
	return res;
};

	//computes all limiters 
	void ComputeLimiters() {
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		limiters.clear();
		for(int i=0; i<cells.size(); i++)
		{
			Cell c = *cells[i];
			limiters[c.GlobalIndex] = ComputeLimiter(c);
		};
		return;
	};

	void Step() {
		//all special calculations
		AditionalStepComputations();

		stepInfo.TimeStep = std::numeric_limits<double>::max();					
		std::vector<Face*> faces = _grid.faces.getLocalNodes();		
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();	
		fluxes.clear();
		#pragma omp for
		for (int i = 0; i<faces.size(); i++) {
			fluxes[faces[i]->GlobalIndex] = std::vector<double>(5,0);
			for (int k = 0; k<5; k++) fluxes[faces[i]->GlobalIndex][k] = 0;		//WID
		};
		vfluxes.clear();
		#pragma omp for
		for (int i = 0; i<faces.size(); i++) {
			vfluxes[faces[i]->GlobalIndex] = std::vector<double>(5,0);
			for (int k = 0; k<5; k++) vfluxes[faces[i]->GlobalIndex][k] = 0;	//WID
		};

		//Time step parameters initialisation
		stepInfo.Residual.resize(5,0);
		//stepInfo.LyambdaC.clear();
		//stepInfo.LyambdaV.clear();
		//stepInfo.LyambdaSum.clear();

		//Compute gradients of conservative variables for each face
		if (IsGradientRequired) 
			ComputeGradients();	
			//ComputeGradientsUniformStruct();

		if(SchemeOrder==2) {
			ComputeConservativeVariablesGradients();
			if(IsLimiterRequired == true) ComputeLimiters();
		};

		//Compute convective fluxes
		ComputeConvectiveFluxes();

		//Compute viscous fluxes		
		if (IsViscousFluxesRequired) 
		{
			ComputeViscousFluxes();			
		};

		//Apply CFL condition
		stepInfo.TimeStep = ComputeMaxTimeStep();	//Blaizek algorithm
		stepInfo.TimeStep *= CFL;

		//Distribute fluxes to cells		
		std::vector<double> sumFlux(5,0);
		//TO DO Remove
		std::vector<double> sumFluxC(5,0);
		std::vector<double> sumFluxV(5,0);

		#pragma omp for
		for (int i = 0; i<cells.size(); i++) {
			Cell c = *cells[i];
			//std::vector<double> sumFluxV = sumFluxCell[c.GlobalIndex];
			bool ExtCell = false;	//TO DO DELETE

			//For each face			
 			for (int k = 0; k<5; k++) sumFlux[k] = 0;
			for (int k = 0; k<5; k++) sumFluxC[k] = 0;
			for (int k = 0; k<5; k++) sumFluxV[k] = 0;

			double P = GetPressure(U[c.GlobalIndex]);	//WID
			
			for (int j = 0; j<c.Faces.size(); j++) {
				Face& f = _grid.faces[c.Faces[j]];
				std::vector<double> flux = fluxes[f.GlobalIndex];// * f.FaceSquare;					
				std::vector<double> vflux = vfluxes[f.GlobalIndex];// * f.FaceSquare;		
				if(f.isExternal) ExtCell = true;
				
				//Consider direction
				//sumFlux += (flux - vflux) * f.FaceSquare;
				if (c.GlobalIndex == f.FaceCell_1) {
					if (!f.isExternal) P = GetPressure(U[f.FaceCell_2]);	//WID
					sumFluxC -= flux * f.FaceSquare;
					sumFluxV += vflux * f.FaceSquare;
					sumFlux -= (flux - vflux) * f.FaceSquare;
				} else {
					if (!f.isExternal) P = GetPressure(U[f.FaceCell_1]);	//WID
					sumFluxC += flux * f.FaceSquare;
					sumFluxV -= vflux * f.FaceSquare;
					sumFlux += (flux - vflux) * f.FaceSquare;
				};
			};
			
			//Add residual
			for (int i = 0; i<5; i++) stepInfo.Residual[i] += sumFlux[i] * sumFlux[i];	
			rouResidual[c.GlobalIndex] = sumFlux[1];
			rovResidual[c.GlobalIndex] = sumFlux[2];
			rowResidual[c.GlobalIndex] = sumFlux[3];
			roResidual[c.GlobalIndex] = sumFlux[0];
			roEResidual[c.GlobalIndex] = sumFlux[4];

			//Add sum flux to cell values
			sumFlux = stepInfo.TimeStep * sumFlux / c.CellVolume;			

			U[c.GlobalIndex].ro += sumFlux[0];
			U[c.GlobalIndex].rou += sumFlux[1];
			U[c.GlobalIndex].rov += sumFlux[2];
			U[c.GlobalIndex].row += sumFlux[3];
			U[c.GlobalIndex].roE += sumFlux[4];

			//uniform external field (RT instability)
			double ro = U[c.GlobalIndex].ro;
			double vx = U[c.GlobalIndex].rou/ro;
			double vy = U[c.GlobalIndex].rov/ro;
			double vz = U[c.GlobalIndex].row/ro;
			U[c.GlobalIndex].rou += g.x*ro*stepInfo.TimeStep;
			U[c.GlobalIndex].rov += g.y*ro*stepInfo.TimeStep;
			U[c.GlobalIndex].row += g.z*ro*stepInfo.TimeStep;
			U[c.GlobalIndex].roE += (g.x*vx + g.y*vy + g.z*vz)*ro*stepInfo.TimeStep;
		};
		
		for (int i = 0; i<5; i++) stepInfo.Residual[i] = sqrt(stepInfo.Residual[i]);	
		totalTime += stepInfo.TimeStep;
	};

	void SetInitialConditionsBlasius(double x0, double U_inf, double rho, double rhoE)
	{
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
	
		for (int i = 0; i<cells.size(); i++) {      
			Cell& cell = *cells[i];
			Vector P = cell.CellCenter;
			ConservativeVariables var;
			//at first density
			var.ro = rho;
			var.roE = rhoE;
			var.row = 0;		//if this component is zero!
			//then velocities
			if(P.x<=x0){
				var.rou = rho*U_inf;
				var.rov = 0;				
			}else{
				double y_BL = P.y*sqrt(U_inf*rho/(medium.Viscosity*(P.x - x0)));
				if(y_BL<3.5)//boundary flow
				{
					var.rou = rho*U_inf*y_BL/3.5;
					var.rov = rho*U_inf*y_BL*P.y/(7*(P.x - x0));
				}else		//outer area
				{
					var.rou = rho*U_inf;
					var.rov = 0;
				};
			};
			U[cell.GlobalIndex] = var;
		};
	
		return;
	};	
	
	//Compute convective fluxes and adjust time step
	void ComputeConvectiveFluxes() {
		std::vector<Face*> faces = _grid.faces.getLocalNodes();		

		//Compute convective flux for each cell face and apply boundary conditions								
		#pragma omp for
		for (int i = 0; i<faces.size(); i++) {
			Face& f = *faces[i];
			if(f.isExternal) fluxes[f.GlobalIndex] = _boundaryConditions[f.BCMarker]->ComputeConvectiveFlux(f);
			else fluxes[f.GlobalIndex] = ComputeConvectiveFlux(f);		
			
			//Additional stability requirement
			double C = 4;
			double lambdaVisc = 0;
			if (IsViscousFluxesRequired) {
				double dX = _grid.cells[f.FaceCell_1].CellVolume / f.FaceSquare;
				lambdaVisc = medium.Viscosity / (dX * U[f.FaceCell_1].ro);
				if (!f.isExternal) {
					double dX2 = _grid.cells[f.FaceCell_2].CellVolume / f.FaceSquare;
					if (lambdaVisc < medium.Viscosity / (dX2 * U[f.FaceCell_2].ro)) lambdaVisc = medium.Viscosity / (dX2 * U[f.FaceCell_2].ro);
				};				
			};

			//Adjust timestep						
			double ts = _grid.cells[f.FaceCell_1].CellVolume/((rSolver.MaxEigenvalue + C * lambdaVisc) * f.FaceSquare);
			if (stepInfo.TimeStep > ts) {
				stepInfo.TimeStep = ts;
			};
			if (!f.isExternal) {				
				ts = _grid.cells[f.FaceCell_2].CellVolume/((rSolver.MaxEigenvalue + C * lambdaVisc) * f.FaceSquare);
				if (stepInfo.TimeStep > ts) stepInfo.TimeStep = ts;
			};

			//Blaizek time step algorithm
			stepInfo.LyambdaC[f.GlobalIndex] = rSolver.MaxEigenvalue*f.FaceSquare;
		};	
	};
	std::vector<double> ComputeConvectiveFlux(Face& f)
	{
		ConservativeVariables UL;
		ConservativeVariables UR;
		switch(SchemeOrder) {
		//first order
		case 1:
			UL = U[f.FaceCell_1];
			if (f.isExternal) {
				UR = GetDummyCellValues(UL, f);
			} else {
				UR = U[f.FaceCell_2];
			};
			break;
		//second order
		case 2:
			Vector dr = f.FaceCenter - _grid.cells[f.FaceCell_1].CellCenter;
			//left reconstraction
			UL = U[f.FaceCell_1];
			std::vector<double> limL = limiters[f.FaceCell_1];	//left limiter values
			UL.ro += (gradCellsRo[f.FaceCell_1]*dr)*limL[0];
			UL.rou += (gradCellsRoU[f.FaceCell_1]*dr)*limL[1];
			UL.rov += (gradCellsRoV[f.FaceCell_1]*dr)*limL[2];
			UL.row += (gradCellsRoW[f.FaceCell_1]*dr)*limL[3];
			UL.roE += (gradCellsRoE[f.FaceCell_1]*dr)*limL[4];
			//right reconstraction
			if (f.isExternal) {
				UR = GetDummyCellValues(UL, f);
			} else {
				dr = f.FaceCenter - _grid.cells[f.FaceCell_2].CellCenter;
				UR = U[f.FaceCell_2];
				std::vector<double> limR = limiters[f.FaceCell_2];	//right limiter values
				UR.ro += (gradCellsRo[f.FaceCell_2]*dr)*limR[0];
				UR.rou += (gradCellsRoU[f.FaceCell_2]*dr)*limR[1];
				UR.rov += (gradCellsRoV[f.FaceCell_2]*dr)*limR[2];
				UR.row += (gradCellsRoW[f.FaceCell_2]*dr)*limR[3];
				UR.roE += (gradCellsRoE[f.FaceCell_2]*dr)*limR[4];
			};
			break;
		};
		//compute appropriate flux
		return rSolver.ComputeFlux( UL,  UR, f);
	};
	double ComputeMaxTimeStep() {
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		double res = std::numeric_limits<double>::max();

		#pragma omp for
		for (int i = 0; i<cells.size(); i++) {
			Cell c = *cells[i];
			double LyambdaSum = 0;
			//special factor
			double C = 4.0/c.CellVolume;
			
			for (int j = 0; j<c.Faces.size(); j++) {
				Face& f = _grid.faces[c.Faces[j]];
				LyambdaSum += stepInfo.LyambdaC[f.GlobalIndex] + C*stepInfo.LyambdaV[f.GlobalIndex];
			};	
			double TS = c.CellVolume/LyambdaSum;
			if(TS<res) res = TS;
		};
		return res;
	};

	//diffusive term functions: turbulent viscosity, heat conductivity and stress tensor
	virtual double GetTurbViscosity(int FaceInd)
	{
		return 0;
	};
	virtual double GetTurbHeatConductivity(int FaceInd)
	{
		return 0;
	};
	virtual Matrix GetAdditionalStresses(int FaceInd)
	{
		Matrix res(3,3);
		return res;
	};

	//Compute viscous flux for each cell face
	void ComputeViscousFluxes() 
	{
		//Obtain all local faces
		std::vector<Face*> faces = _grid.faces.getLocalNodes();				
		std::vector<double> sumFlux(5,0);

		//Compute viscous flux for each face
		#pragma omp for
		for (int i = 0; i<faces.size(); i++) {
			Face& f = *faces[i];		

			std::vector<double> vflux(5, 0);
			for (int k = 0; k<5; k++) vflux[k] = 0;
			
			Vector CL = _grid.cells[f.FaceCell_1].CellCenter;
			Vector CR;
			ConservativeVariables UL = U[f.FaceCell_1]; 
            ConservativeVariables UR; 
            if (f.isExternal) { 
                UR = GetDummyCellValues(UL, f); 
				CR = 2 * (f.FaceCenter - _grid.cells[f.FaceCell_1].CellCenter) + _grid.cells[f.FaceCell_1].CellCenter;
            } else { 
                UR = U[f.FaceCell_2]; 
				CR = _grid.cells[f.FaceCell_2].CellCenter;
            };
			
			//Here we compute actual viscocity thermal conductivity and stresses
	
			//get turbulent viscosity
			double viscosityTurb = GetTurbViscosity(f.GlobalIndex);
			double KTurb = GetTurbHeatConductivity(f.GlobalIndex);
			Matrix Stresses = GetAdditionalStresses(f.GlobalIndex);
			
			//double viscosity = medium.Viscosity;	
			double viscosity = medium.Viscosity;
			double s_viscosity = -2.0/3.0 * viscosity;		//second viscosity STOKES LAW
			viscosity += viscosityTurb;			//for bussinesk hypothesis
	
			//required variables on face	
			ConservativeVariables UAvg;
			//for (int k = 0; k<5; k++) UAvg[k] = 0.5*(UL[k] + UR[k]);
			double dL = abs(f.FaceNormal * (CL - f.FaceCenter));		//WID
			double dR = abs(f.FaceNormal * (CR - f.FaceCenter));				
			for (int k = 0; k<5; k++) UAvg[k] = (UL[k] * dR + UR[k] * dL) / (dL + dR);		
			double rho = GetDensity(UAvg);
			double u = GetVelocityX(UAvg);
			double v = GetVelocityY(UAvg);
			double w = GetVelocityZ(UAvg);			

			//gradients of required variables
			Vector du, dv, dw, dT;

			//New try
			du = gradFacesU[f.GlobalIndex];
			dv = gradFacesV[f.GlobalIndex];
			dw = gradFacesW[f.GlobalIndex];
			dT = gradFacesT[f.GlobalIndex];			

			double R = (medium.Gamma - 1.0) * medium.Cv;		//specific gas constant
			double K = medium.ThermalConductivity + KTurb;		//effective heat conduction coefficient										

			//stress elements
			double tau_diagonal = s_viscosity*(du.x + dv.y + dw.z);
			double tau_xx = tau_diagonal + 2*viscosity*du.x + Stresses[0][0];
			double tau_yy = tau_diagonal + 2*viscosity*dv.y + Stresses[1][1];
			double tau_zz = tau_diagonal + 2*viscosity*dw.z + Stresses[2][2];
			double tau_xy = viscosity*(du.y + dv.x) + Stresses[0][1];
			double tau_xz = viscosity*(du.z + dw.x) + Stresses[0][2];
			double tau_yz = viscosity*(dv.z + dw.y) + Stresses[1][2];						

			//work of viscous stresses and heat conduction
			Vector Thetta;
			Thetta.x = u*tau_xx + v*tau_xy + w*tau_xz + K*dT.x;
			Thetta.y = u*tau_xy + v*tau_yy + w*tau_yz + K*dT.y;
			Thetta.z = u*tau_xz + v*tau_yz + w*tau_zz + K*dT.z;		

			//viscous fluxes
 			vflux[1] = f.FaceNormal.x*tau_xx + f.FaceNormal.y*tau_xy + f.FaceNormal.z*tau_xz;			
			vflux[2] = f.FaceNormal.x*tau_xy + f.FaceNormal.y*tau_yy + f.FaceNormal.z*tau_yz;			
			vflux[3] = f.FaceNormal.x*tau_xz + f.FaceNormal.y*tau_yz + f.FaceNormal.z*tau_zz;			
			vflux[4] = f.FaceNormal*Thetta;				

			//Add viscous fluxes		
			sumFlux = abs(vflux);
			std::vector<double> cflux = fluxes[f.GlobalIndex];

			//TO DO remove			
			vfluxes[f.GlobalIndex] = vflux;

			//Blaizek time step algorithm
			stepInfo.LyambdaV[f.GlobalIndex] = max(4.0/3.0, medium.Gamma)*(medium.ThermalConductivity/(medium.Cv*medium.Gamma*rho))*f.FaceSquare*f.FaceSquare;
		};
		
	};	

	//Function to compute gradients at every cell
	void ComputeGradients() {	
		//Clear all previous data TO DO REMOVE
		gradCellsU.clear();
		gradCellsV.clear();
		gradCellsW.clear();
		gradCellsT.clear();
		gradFacesU.clear();
		gradFacesV.clear();
		gradFacesW.clear();
		gradFacesT.clear();						
		
		//Compute required gradients in cells
		ComputeFunctionGradient(gradCellsU, U, &Model<RiemannSolver>::GetVelocityX);
		ComputeFunctionGradient(gradCellsV, U, &Model<RiemannSolver>::GetVelocityY);
		ComputeFunctionGradient(gradCellsW, U, &Model<RiemannSolver>::GetVelocityZ);
		ComputeFunctionGradient(gradCellsT, U, &Model<RiemannSolver>::GetTemperature);

		//Compute gradients at face centers
		//Obtain all local faces
		std::vector<Face*> faces = _grid.faces.getLocalNodes();		
		#pragma omp for
		for (int i = 0; i<faces.size(); i++) {
			Face& face = *faces[i];
			
			//Preparation
			Vector gU;
			Vector gV;
			Vector gW;
			Vector gT;
			Vector CL = _grid.cells[face.FaceCell_1].CellCenter;
			Vector CR;
			ConservativeVariables UL = U[face.FaceCell_1]; 
            ConservativeVariables UR; 
            if (face.isExternal) { 
                UR = GetDummyCellValues(UL, face); 
				CR = 2 * (face.FaceCenter - _grid.cells[face.FaceCell_1].CellCenter) + _grid.cells[face.FaceCell_1].CellCenter;
            } else { 
                UR = U[face.FaceCell_2]; 
				CR = _grid.cells[face.FaceCell_2].CellCenter;
            };

			//Compute average of gradient on face
			if (face.isExternal) { //If external just take simple gradient from nearest cell
				gU = gradCellsU[face.FaceCell_1];
				gV = gradCellsV[face.FaceCell_1];
				gW = gradCellsW[face.FaceCell_1];
				gT = gradCellsT[face.FaceCell_1];
				bool isWall = _wallBoundaryMarkers.find(face.BCMarker) != _wallBoundaryMarkers.end();
				// TO DO REMOVE KOSTIL WID
				if (isWall) { //Wall additional correction
					gU.x = 0;
					gU.y = U[face.FaceCell_1].rou / (U[face.FaceCell_1].ro*_grid.cells[face.FaceCell_1].CellCenter.y);
					gV.x = 0;
					gV.y = U[face.FaceCell_1].rov / (U[face.FaceCell_1].ro*_grid.cells[face.FaceCell_1].CellCenter.y);
					gT.y = 0;
				};
				//
				

			} else { //Linear averaging						//WID я тут исправил интерполяцию		
				double dL = abs(face.FaceNormal * (CL - face.FaceCenter));
				double dR = abs(face.FaceNormal * (CR - face.FaceCenter));				
				gU = (gradCellsU[face.FaceCell_1] * dR + gradCellsU[face.FaceCell_2] * dL) / (dL + dR);				
				gV = (gradCellsV[face.FaceCell_1] * dR + gradCellsV[face.FaceCell_2] * dL) / (dL + dR);				
				gW = (gradCellsW[face.FaceCell_1] * dR + gradCellsW[face.FaceCell_2] * dL) / (dL + dR);				
				gT = (gradCellsT[face.FaceCell_1] * dR + gradCellsT[face.FaceCell_2] * dL) / (dL + dR);				
			};

			//Gradient correction
			Vector dVelocity = GetVelocity(UR) - GetVelocity(UL);
			double dT = GetTemperature(UR) - GetTemperature(UL);
			Vector dr = CR - CL;
			Vector t = dr / dr.mod();
			double dl = dr.mod();
			gU = gU - (t * gU - (dVelocity.x)/(dl)) * t;			
			gV = gV - (t * gV - (dVelocity.y)/(dl)) * t;			
			gW = gW - (t * gW - (dVelocity.z)/(dl)) * t;	
			gT = gT - (t * gT - (dT)/(dl)) * t;	

			//Eventually
			gradFacesU[face.GlobalIndex] = gU;
			gradFacesV[face.GlobalIndex] = gV;
			gradFacesW[face.GlobalIndex] = gW;
			gradFacesT[face.GlobalIndex] = gT;
		};
		
	};
	//Function to compute gradients at every cell
	void ComputeGradientsUniformStruct() {	
		//Clear all previous data
		gradCellsU.clear();
		gradCellsV.clear();
		gradCellsW.clear();
		gradCellsT.clear();
		gradFacesU.clear();
		gradFacesV.clear();
		gradFacesW.clear();
		gradFacesT.clear();

		//compute gradients over all cells
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		for(int i=0; i<cells.size(); i++)
		{
			Cell c = *cells[i];
			ConservativeVariables U_bottom, U_top, U_left, U_right;
			Vector R_bottom, R_top, R_left, R_right;
			for(int j=0; j<c.Faces.size();j++)
			{
				ConservativeVariables Ub;
				Face f = _grid.faces[c.Faces[j]];
				if(f.isExternal) Ub = GetDummyCellValues(U[c.GlobalIndex], f);
				else {
					if(c.GlobalIndex==f.FaceCell_1) Ub = U[f.FaceCell_2];
					else Ub = U[f.FaceCell_1];
				};
				//define relative position
				if(abs(f.FaceNormal.x)>abs(f.FaceNormal.y)) {
					if(c.CellCenter.x>f.FaceCenter.x) {
						//left
						U_left = Ub;
						R_left = f.FaceCenter;
					}else{
						//right
						U_right = Ub;
						R_right = f.FaceCenter;
					};
				}else {
					if(c.CellCenter.y>f.FaceCenter.y) {
						//bottom
						U_bottom = Ub;
						R_bottom = f.FaceCenter;
					}else{
						//top
						U_top = Ub;
						R_top = f.FaceCenter;
					};
				};
			};
			//compute gradients in cell c
			//X velocity
			gradCellsU[c.GlobalIndex].z = 0;
			gradCellsU[c.GlobalIndex].x = 0.5*(GetVelocityX(U_right) - GetVelocityX(U_left))/(R_right - R_left).x;
			gradCellsU[c.GlobalIndex].y = 0.5*(GetVelocityX(U_top) - GetVelocityX(U_bottom))/(R_top - R_bottom).y;
			//Y velocity
			gradCellsV[c.GlobalIndex].z = 0;
			gradCellsV[c.GlobalIndex].x = 0.5*(GetVelocityY(U_right) - GetVelocityY(U_left))/(R_right - R_left).x;
			gradCellsV[c.GlobalIndex].y = 0.5*(GetVelocityY(U_top) - GetVelocityY(U_bottom))/(R_top - R_bottom).y;
			//Z velocity
			gradCellsW[c.GlobalIndex].z = 0;
			gradCellsW[c.GlobalIndex].x = 0;
			gradCellsW[c.GlobalIndex].y = 0;
			//Temperature
			gradCellsT[c.GlobalIndex].z = 0;
			gradCellsT[c.GlobalIndex].x = 0.5*(GetTemperature(U_right) - GetTemperature(U_left))/(R_right - R_left).x;
			gradCellsT[c.GlobalIndex].y = 0.5*(GetTemperature(U_top) - GetTemperature(U_bottom))/(R_top - R_bottom).y;
		};
		
		//Compute gradients at face centers
		//Obtain all local faces
		std::vector<Face*> faces = _grid.faces.getLocalNodes();		
		#pragma omp for
		for (int i = 0; i<faces.size(); i++) {
			Face& face = *faces[i];
			
			//Preparation
			Vector gU;
			Vector gV;
			Vector gW;
			Vector gT;
			ConservativeVariables UL = U[face.FaceCell_1];
			Vector CL = _grid.cells[face.FaceCell_1].CellCenter;
            ConservativeVariables UR;
			Vector CR;
            if (face.isExternal) { 
                UR = GetDummyCellValues(UL, face);
				CR = 2*(face.FaceCenter - CL) + CL;
			} else { 
                UR = U[face.FaceCell_2];
				CR = _grid.cells[face.FaceCell_2].CellCenter;
            };
			//compute gradient at faces
			if(abs(face.FaceNormal.x)>abs(face.FaceNormal.y))	//vertical face
			{
				gU.x = (GetVelocityX(UR) - GetVelocityX(UL))/(CR - CL).x;
				gV.x = (GetVelocityY(UR) - GetVelocityY(UL))/(CR - CL).x;
				gW.x = 0;
				gT.x = (GetTemperature(UR) - GetTemperature(UL))/(CR - CL).x;

				if(face.isExternal) {
					gU.y = gradCellsU[face.FaceCell_1].y;
					gV.y = gradCellsV[face.FaceCell_1].y;
					gW.y = 0;
					gT.y = gradCellsT[face.FaceCell_1].y;
				}else {
					gU.y = 0.5*(gradCellsU[face.FaceCell_1].y + gradCellsU[face.FaceCell_2].y);
					gV.y = 0.5*(gradCellsV[face.FaceCell_1].y + gradCellsV[face.FaceCell_2].y);
					gW.y = 0;
					gT.y = 0.5*(gradCellsT[face.FaceCell_1].y + gradCellsT[face.FaceCell_2].y);
				};
			} else //horizontal face
			{	
				gU.y = (GetVelocityX(UR) - GetVelocityX(UL))/(CR - CL).y;
				gV.y = (GetVelocityY(UR) - GetVelocityY(UL))/(CR - CL).y;
				gW.y = 0;
				gT.y = (GetTemperature(UR) - GetTemperature(UL))/(CR - CL).y;

				if(face.isExternal) {
					gU.x = gradCellsU[face.FaceCell_1].y;
					gV.x = gradCellsV[face.FaceCell_1].y;
					gW.x = 0;
					gT.x = gradCellsT[face.FaceCell_1].y;
				} else {
					gU.x = 0.5*(gradCellsU[face.FaceCell_1].x + gradCellsU[face.FaceCell_2].x);
					gV.x = 0.5*(gradCellsV[face.FaceCell_1].x + gradCellsV[face.FaceCell_2].x);
					gW.x = 0;
					gT.x = 0.5*(gradCellsT[face.FaceCell_1].x + gradCellsT[face.FaceCell_2].x);
				};
			};
			//Z components
			gU.z = 0;
			gV.z = 0;
			gW.z = 0;
			gT.z = 0;

			//Eventually
			gradFacesU[face.GlobalIndex] = gU;
			gradFacesV[face.GlobalIndex] = gV;
			gradFacesW[face.GlobalIndex] = gW;
			gradFacesT[face.GlobalIndex] = gT;
		};
		
	};
	//compute conservative veriables gradients in all cells
	void ComputeConservativeVariablesGradients()
	{
		ComputeFunctionGradient(gradCellsRo, U, &Model<RiemannSolver>::GetDensity);
		ComputeFunctionGradient(gradCellsRoU, U, &Model<RiemannSolver>::GetRoU);
		ComputeFunctionGradient(gradCellsRoV, U, &Model<RiemannSolver>::GetRoV);
		ComputeFunctionGradient(gradCellsRoW, U, &Model<RiemannSolver>::GetRoW);
		ComputeFunctionGradient(gradCellsRoE, U, &Model<RiemannSolver>::GetRoE);
	};

	//Compute wall distances for every cell
	void ComputeWallDistances() {		
		//If no walls defined
		if (_wallFaces.size() == 0) {
			std::cout<<"No walls defined. Distance computation is omitted.\n";
			return;
		};		

		// Build alglib kdtree for all wall face centers
		alglib::ae_int_t nx = 3;
		alglib::ae_int_t ny = 1;
		alglib::ae_int_t normtype = 2;		
		alglib::real_2d_array a;  		
		a.setlength(_wallFaces.size(), nx + ny);		

		int i = 0;
		for (std::map<int, FaceWallInfo>::iterator it = _wallFaces.begin(); it != _wallFaces.end(); it++) {			
			// Clear all wall info for faces
			it->second.layerCells.clear();

			//Fill in kdtree information
			a[i][0] = _grid.faces[it->first].FaceCenter.x;
			a[i][1] = _grid.faces[it->first].FaceCenter.y;
			a[i][2] = _grid.faces[it->first].FaceCenter.z;
			a[i][3] = _grid.faces[it->first].GlobalIndex;
			i++;
		};
	
		// Build tree	
		alglib::kdtreebuild(a, nx, ny, normtype, kdtWall);		

		//Now for each local cell obtain distance to closest wall face and other properties
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		alglib::real_1d_array x; 
		x.setlength(nx);
		alglib::real_2d_array xy; 		
		xy.setlength(1, nx+ny);
		for (int i = 0; i < cells.size(); i++) {			
			x[0] = cells[i]->CellCenter.x;
			x[1] = cells[i]->CellCenter.y;
			x[2] = cells[i]->CellCenter.z;			
			alglib::ae_int_t isFound = alglib::kdtreequeryknn(kdtWall, x, 1);
			alglib::kdtreequeryresultsxy(kdtWall, xy);						
			Vector radiusVector = Vector(xy[0][0], xy[0][1], xy[0][2]);			
			radiusVector -= cells[i]->CellCenter;
			int faceIndex = xy[0][3];
			_wallInfo[cells[i]->GlobalIndex] = CellWallInfo();		
			_wallInfo[cells[i]->GlobalIndex].wallFaceIndex = faceIndex;
			_wallInfo[cells[i]->GlobalIndex].distance = radiusVector.mod();
			_wallFaces[faceIndex].layerCells.push_back(cells[i]->GlobalIndex);
		};

		//compute angles
		for (int i = 0; i < cells.size(); i++)
		{
			Face f = _grid.faces[_wallInfo[cells[i]->GlobalIndex].wallFaceIndex];
			Vector FaceCell = cells[i]->CellCenter - f.FaceCenter;
			double Cos = abs(FaceCell*f.FaceNormal)/(f.FaceNormal.mod()*FaceCell.mod());
			_wallInfo[cells[i]->GlobalIndex].angle  = acos(Cos);
		};

		//For each wall face sort cells ascendig by distance
		DistanceSorting();
	};

	//Compute hight of boundary layer
	void ComputeBoundaryLayerHeight(double angle)
	{
		for (std::map<int, FaceWallInfo>::iterator it = _wallFaces.begin(); it != _wallFaces.end(); it++)
		{
			int FaceInd = it->first;	//index of during wall face
			//list of cells up to our face
			std::vector<int> CellSequence(0);

			//Get appropriate cells by angle and fill in the list
			for(int i=0; i<_wallFaces[FaceInd].layerCells.size(); i++) {
			CellWallInfo& wInfo = _wallInfo[_wallFaces[FaceInd].layerCells[i]];
			if(wInfo.angle<angle) CellSequence.push_back(_wallFaces[FaceInd].layerCells[i]);
			};

			//get first local maximum of appropriate velocity component as the top of boundary layer
			int UpCell;
			double U_tan_max = 0;
			for(int i=0; i<CellSequence.size(); i++)
			{
				Vector Velocity(U[CellSequence[i]].rou/U[CellSequence[i]].ro, U[CellSequence[i]].rov/U[CellSequence[i]].ro, U[CellSequence[i]].row/U[CellSequence[i]].ro);
				double U_tan = (Velocity - (Velocity*_grid.faces[FaceInd].FaceNormal)*_grid.faces[FaceInd].FaceNormal).mod();	//tangential velocity

				if(U_tan > U_tan_max)
				{
					U_tan_max = U_tan;	
					UpCell = CellSequence[i];
				}else break;
			};

			//write result
			_boundaryLayerHeight[FaceInd] = _wallInfo[UpCell].distance;
		};
	};

	//object for compare 2 cells for quick sort	
	struct ComparerCellWallDist {
	private:
		std::map<int, CellWallInfo>& wallInfo;
	public:
		ComparerCellWallDist(std::map<int, CellWallInfo>& _wI) : wallInfo(_wI) {};
		bool operator() (int i,int j) {
			return wallInfo[i].distance < wallInfo[j].distance;			
		};		
	};

	//Sort cells layer from wall
	void DistanceSorting()
	{
		//initialize comparer object
		ComparerCellWallDist cmp(_wallInfo);

		//loop on all wall faces
		for (std::map<int, FaceWallInfo>::iterator it = _wallFaces.begin(); it != _wallFaces.end(); it++)
		{			
			std::sort(it->second.layerCells.begin(), it->second.layerCells.end(), cmp);
			for (int j = 1; j<it->second.layerCells.size(); j++) {
				if (_wallInfo[it->second.layerCells[j-1]].distance > _wallInfo[it->second.layerCells[j]].distance) std::cout<<"Sort fail\n";	//WID - лишняя проверка
			};
		};	
	};
	
	//compute displacement thickness with cells soarting by angle
	double ComputeDisplacementThickness(int wallFaceInd, double U_inf)
	{
		//layer cells must be sorted in ascending order by distance in some angular range
		std::vector<int> CellSequence(0);
		for(int i=0; i<_wallFaces[wallFaceInd].layerCells.size(); i++) {
			CellWallInfo& wInfo = _wallInfo[_wallFaces[wallFaceInd].layerCells[i]];
			if(wInfo.angle<SoartingAngle) {
				CellSequence.push_back(_wallFaces[wallFaceInd].layerCells[i]);
			};
		};
		
		//write required function
		function f(CellSequense.size()+1);
		f.x[0] = 0;		//values on wall
		f.value[0] = 1.0;
		for(int i=0; i<CellSequense.size(); i++)
		{
			f.x[i+1] = _wallInfo[CellSequence[i]].distance;
			f.value[i+1] = 1 - U[CellSequence[i]].rou/(U[CellSequence[i]].ro*U_inf);
		};
		//compute integral by trapezium formula
		return f.TrapeziumIntegral();
	};
	//compute displacement thickness by velocity distribution from wall to top of boundary layer
	double ComputeDisplacementThickness(const function& velosity)
	{
		//write required function
		function f = velosity;
		double U_top = velosity.value[velosity.size - 1];
		for(int i=0; i<f.size; i++)
			f.value[i] = 1.0 - velosity.value[i]/U_top;

		//compute integral by trapezium formula
		return f.TrapeziumIntegral();
	};

	//Boundary conditions main processing function
	ConservativeVariables GetDummyCellValues(ConservativeVariables U, const Face& face) {
		//Determine boundary condition type
		BoundaryCondition* bc = _boundaryConditions[face.BCMarker];
		return bc->getDummyValues(U, face); 
	};	

	//Conversion functions and postprocessing
	ConservativeVariables PrimitiveToConservativeVariables(Vector v, double p, double T, MediumProperties medium) {
		ConservativeVariables res;
		double e = medium.Cv * T;
		double v2 = v*v; 
		res.ro = p / (e * (medium.Gamma - 1));
		res.rou = v.x * res.ro;
		res.rov = v.y * res.ro;
		res.row = v.z * res.ro;
		res.roE = (e + v2 / 2.0)* res.ro;
		return res;
	};	
	Vector GetVelocity(const ConservativeVariables& celldata) {		
		double ro = celldata.ro;
		double vx = celldata.rou/celldata.ro;
		double vy = celldata.rov/celldata.ro;
		double vz = celldata.row/celldata.ro;
		return Vector(vx, vy, vz);
	};
	double GetVelocityX(const ConservativeVariables& celldata) {		
		double ro = celldata.ro;
		double vx = celldata.rou/celldata.ro;		
		return vx;
	};
	double GetVelocityY(const ConservativeVariables& celldata) {		
		double ro = celldata.ro;
		double vy = celldata.rov/celldata.ro;		
		return vy;
	};
	double GetVelocityZ(const ConservativeVariables& celldata) {		
		double ro = celldata.ro;
		double vz = celldata.row/celldata.ro;		
		return vz;
	};
	double GetPressure(const ConservativeVariables& celldata) {
		double ro = celldata.ro;
		double vx = celldata.rou/celldata.ro;
		double vy = celldata.rov/celldata.ro;
		double vz = celldata.row/celldata.ro;
		double E = celldata.roE/celldata.ro;
		double P = (medium.Gamma-1.0) * ro * (E - (vx*vx+vy*vy+vz*vz)/2.0);
		return P;
	};
	double GetTemperature(const ConservativeVariables& celldata) {
		double ro = celldata.ro;
		double vx = celldata.rou/celldata.ro;
		double vy = celldata.rov/celldata.ro;
		double vz = celldata.row/celldata.ro;
		double E = celldata.roE/celldata.ro;
		double T = (E - (vx*vx+vy*vy+vz*vz)/2.0) / (medium.Cv);
		return T;
	};
	double GetSoundSpeed(const ConservativeVariables& celldata) {
		double ro = celldata.ro;
		double vx = celldata.rou/celldata.ro;
		double vy = celldata.rov/celldata.ro;
		double vz = celldata.row/celldata.ro;
		double E = celldata.roE/celldata.ro;
		double P = (medium.Gamma-1.0) * ro * (E - (vx*vx+vy*vy+vz*vz)/2.0);
		double c = sqrt(medium.Gamma * P / ro);
		return c;
	};
	//conservative variables
	double GetDensity(const ConservativeVariables& celldata) {		
		double ro = celldata.ro;		
		return ro;
	};
	double GetRoU(const ConservativeVariables& celldata) {    
	 return celldata.rou;
	};
	double GetRoV(const ConservativeVariables& celldata) {    
	return celldata.rov;
	};
	double GetRoW(const ConservativeVariables& celldata) {    
	 return celldata.row;
	};
	double GetRoE(const ConservativeVariables& celldata) {    
	 return celldata.roE;
	};

	void SaveWallPressureDistribution(std::string fname) {
		std::ofstream ofs(fname);
		ofs<<std::scientific;

		//Header
		ofs<<"VARIABLES= \"X\"\n";
		ofs<<"\"Y\"\n";
		ofs<<"\"Z\"\n";
		ofs<<"\"P\"\n";
		ofs<<"ZONE\n";

		for (std::map<int, FaceWallInfo>::iterator it = _wallFaces.begin(); it != _wallFaces.end(); it++) {
			FaceWallInfo& wInfo = it->second;
			Face& face = _grid.faces[wInfo.wallFaceIndex];
			double P = GetPressure(U[face.FaceCell_1]);
			ofs<<face.FaceCenter.x<<" ";
			ofs<<face.FaceCenter.y<<" ";
			ofs<<face.FaceCenter.z<<" ";
			ofs<<P<<"\n";
		};			
		
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		std::vector<ConservativeVariables*> data = U.getLocalNodes();

		ofs.close();
	};

	//WID
	virtual void SaveToTecplot25D(std::string fname) {
		std::ofstream ofs(fname);
		ofs<<std::scientific;

		//Access local nodes, faces, cells and flow data
		std::vector<Node*> nodesOriginal = _grid.nodes.getLocalNodes();
		std::vector<Face*> facesOriginal = _grid.faces.getLocalNodes();
		std::vector<Cell*> cellsOriginal = _grid.cells.getLocalNodes();
		std::vector<ConservativeVariables*> data = U.getLocalNodes();

		//Collect everything from Z=0 plane
		std::map<int,int> toNaturalIndex;			
		int counter = 1;
		std::vector<Node*> nodes;
		std::vector<Cell*> cells;
		for (int i = 0; i<nodesOriginal.size(); i++) {
			if (nodesOriginal[i]->P.z < 1e-2) {
				nodes.push_back(nodesOriginal[i]);
				//Map all node global indexes to natural numbers
				toNaturalIndex[nodesOriginal[i]->GlobalIndex] = counter++;
			};
		};

		std::map<int, std::vector<int>> cellPlaneFaceNodes;
		for (int i = 0; i<cellsOriginal.size(); i++) {
			cells.push_back(cellsOriginal[i]);
			//Find face in Z=0 plane
			bool isFound = false;
			Cell& c = *cellsOriginal[i];
			for (int j = 0; j<c.Faces.size(); j++) {
				if (_grid.faces[c.Faces[j]].FaceCenter.z < 1e-2) {
					cellPlaneFaceNodes[c.GlobalIndex] = _grid.faces[c.Faces[j]].FaceNodes;
					isFound = true;
					break;
				};				
			};
			if (!isFound) {
				std::cout<<"Face in plane Z=0 not found\n";
			};
			//if (toNaturalIndex.find(cellPlaneFaceNodes[cellsOriginal[i]->GlobalIndex][0]) == toNaturalIndex.end()) std::cout<<"Node not from Z=0 plane\n";
			//if (toNaturalIndex.find(cellPlaneFaceNodes[cellsOriginal[i]->GlobalIndex][1]) == toNaturalIndex.end()) std::cout<<"Node not from Z=0 plane\n";
			//if (toNaturalIndex.find(cellPlaneFaceNodes[cellsOriginal[i]->GlobalIndex][2]) == toNaturalIndex.end()) std::cout<<"Node not from Z=0 plane\n";
			//if (toNaturalIndex.find(cellPlaneFaceNodes[cellsOriginal[i]->GlobalIndex][3]) == toNaturalIndex.end()) std::cout<<"Node not from Z=0 plane\n";
		};

		//TO DO unify		
		//Header
		ofs<<"VARIABLES= \"X\", \"Y\", \"Rho\", \"u\", \"v\", \"w\", \"E\", \"T\", \"P\"";
		ofs<<", \"PStagnation\"";
		ofs<<", \"distance\"";
		ofs<<"\n";
			
		ofs<<"ZONE T=\"D\"\n";
		ofs<<"N=" << nodes.size() << ", E=" << cells.size() <<", F=FEBLOCK, ET=QUADRILATERAL\n";

		ofs<<"VARLOCATION = (NODAL, NODAL, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED";
		ofs<<", CELLCENTERED";
		ofs<<", CELLCENTERED";
		ofs<<")\n";

		//Nodes coordinates
		//X
		for (int i = 0; i<nodes.size(); i++) {
			ofs<<nodes[i]->P.x<<"\n";
		};

		//Y
		for (int i = 0; i<nodes.size(); i++) {
			ofs<<nodes[i]->P.y<<"\n";
		};

		////Solution data
		////Rho
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<U[idx].ro<<"\n";
		};
				
		//u
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<U[idx].rou/U[idx].ro<<"\n";
		};
		//v
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<U[idx].rov/U[idx].ro<<"\n";		
		};
		//w
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<U[idx].row/U[idx].ro<<"\n";
		};
		//E
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<U[idx].roE/U[idx].ro<<"\n";
		};
		//T
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			double ro = U[idx].ro;
			double vx = U[idx].rou/U[idx].ro;
			double vy = U[idx].rov/U[idx].ro;
			double vz = U[idx].row/U[idx].ro;
			double E = U[idx].roE/U[idx].ro;
			double T = (E - (vx*vx+vy*vy+vz*vz)/2.0) / (medium.Cv);
			ofs<<T<<"\n";
		};
		//P
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			double ro = U[idx].ro;
			double vx = U[idx].rou/U[idx].ro;
			double vy = U[idx].rov/U[idx].ro;
			double vz = U[idx].row/U[idx].ro;
			double E = U[idx].roE/U[idx].ro;
			double P = (medium.Gamma-1.0) * ro * (E - (vx*vx+vy*vy+vz*vz)/2.0);
			ofs<<P<<"\n";
		};	

		//PStagnation
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			double ro = U[idx].ro;
			double vx = U[idx].rou/U[idx].ro;
			double vy = U[idx].rov/U[idx].ro;
			double vz = U[idx].row/U[idx].ro;
			double E = U[idx].roE/U[idx].ro;
			double PStagnation = (medium.Gamma-1.0) * U[idx].roE;
			ofs<<PStagnation<<"\n";
		};	

		//distance
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			if (_wallInfo.size() != 0) {
				ofs<<_wallInfo[idx].distance<<"\n";
			} else {
				ofs<<0<<"\n";
			};
		};

		//Connectivity list for each cell
		for (int i = 0; i<cells.size(); i++) {
			ofs<<toNaturalIndex[cellPlaneFaceNodes[cells[i]->GlobalIndex][0]]<<" "
				<<toNaturalIndex[cellPlaneFaceNodes[cells[i]->GlobalIndex][1]]<<" "
				<<toNaturalIndex[cellPlaneFaceNodes[cells[i]->GlobalIndex][2]]<<" "
				<<toNaturalIndex[cellPlaneFaceNodes[cells[i]->GlobalIndex][3]]<<"\n";
		};		

		ofs.close();
		return;
	};

	//Save current model state to techplot file
	virtual void SaveToTechPlotWithGrad(std::string fname) {
		std::ofstream ofs(fname);
		ofs<<std::scientific;
		if (_grid.gridInfo.GridDimensions == 1) {
			//Header
			ofs<<"VARIABLES= \"x\"\n";
			ofs<<"\"ro\"\n";
			ofs<<"\"u\"\n";
			ofs<<"\"P\"\n";
			ofs<<"\"T\"\n";
			ofs<<"ZONE\n";

			std::vector<Cell*> cells = _grid.cells.getLocalNodes();
			for (int i = 0; i<cells.size(); i++) {			
				Cell& c = *cells[i];								
				ofs<<c.CellCenter.x<<"\n";
				ofs<<U[c.GlobalIndex].ro<<"\n";
				ofs<<GetVelocity(U[c.GlobalIndex]).x<<"\n";
				ofs<<GetPressure(U[c.GlobalIndex])<<"\n";
				ofs<<GetTemperature(U[c.GlobalIndex])<<"\n";
				ofs<<"\n";
			};		
		};

		if (_grid.gridInfo.GridDimensions == 2) {
			//TO DO unify		
			//Header

			//Access local nodes, faces, cells and flow data
			std::vector<Node*> nodes = _grid.nodes.getLocalNodes();
			std::vector<Face*> faces = _grid.faces.getLocalNodes();
			std::vector<Cell*> cells = _grid.cells.getLocalNodes();
			std::vector<ConservativeVariables*> data = U.getLocalNodes();

			ofs<<"VARIABLES= \"X\", \"Y\", \"Rho\", \"u\", \"v\", \"w\", \"E\", \"T\", \"P\"";
			ofs<<", \"PStagnation\"";
			ofs<<", \"distance\"";
			ofs<<", \"rouResidual\"";
			ofs<<", \"rovResidual\"";
			ofs<<", \"rowResidual\"";
			ofs<<", \"roResidual\"";
			ofs<<", \"roEResidual\"";
			ofs<<", \"roU\"";
			ofs<<", \"roV\"";
			ofs<<", \"dU_dx\"";
			ofs<<", \"dU_dy\"";
			ofs<<"\n";
			
			ofs<<"ZONE T=\"D\"\n";
			ofs<<"N=" << _grid.nodes.size() << ", E=" << _grid.cells.size() <<", F=FEBLOCK, ";
			if (cells[0]->CGNSType == QUAD_4) {
				ofs<<"ET=QUADRILATERAL\n";
			};
			if (cells[0]->CGNSType == TRI_3) {
				ofs<<"ET=TRIANGLE\n";
			};

			ofs<<"VARLOCATION = (NODAL, NODAL, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<")\n";

			//Map all node global indexes to natural numbers
			std::map<int,int> toNaturalIndex;
			std::set<int> nodeIndexes = _grid.nodes.getAllIndexes();
			int counter = 1;
			for (std::set<int>::iterator it = nodeIndexes.begin(); it != nodeIndexes.end(); it++) toNaturalIndex[*it] = counter++;			

			//Nodes coordinates
			//X
			for (int i = 0; i<nodes.size(); i++) {
				ofs<<nodes[i]->P.x<<"\n";
			};

			//Y
			for (int i = 0; i<nodes.size(); i++) {
				ofs<<nodes[i]->P.y<<"\n";
			};

			////Solution data
			////Rho
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].ro<<"\n";
			};
				
			//u
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rou/U[idx].ro<<"\n";
			};
			//v
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rov/U[idx].ro<<"\n";		
			};
			//w
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].row/U[idx].ro<<"\n";
			};
			//E
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].roE/U[idx].ro<<"\n";
			};
			//T
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				double ro = U[idx].ro;
				double vx = U[idx].rou/U[idx].ro;
				double vy = U[idx].rov/U[idx].ro;
				double vz = U[idx].row/U[idx].ro;
				double E = U[idx].roE/U[idx].ro;
				double T = (E - (vx*vx+vy*vy+vz*vz)/2.0) / (medium.Cv);
				ofs<<T<<"\n";
			};
			//P
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				double ro = U[idx].ro;
				double vx = U[idx].rou/U[idx].ro;
				double vy = U[idx].rov/U[idx].ro;
				double vz = U[idx].row/U[idx].ro;
				double E = U[idx].roE/U[idx].ro;
				double P = (medium.Gamma-1.0) * ro * (E - (vx*vx+vy*vy+vz*vz)/2.0);
				ofs<<P<<"\n";
			};	

			//PStagnation
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				double ro = U[idx].ro;
				double vx = U[idx].rou/U[idx].ro;
				double vy = U[idx].rov/U[idx].ro;
				double vz = U[idx].row/U[idx].ro;
				double E = U[idx].roE/U[idx].ro;
				double PStagnation = (medium.Gamma-1.0) * U[idx].roE;
				ofs<<PStagnation<<"\n";
			};	

			//distance
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				if (_wallInfo.size() != 0) {
					ofs<<_wallInfo[idx].distance<<"\n";
				} else {
					ofs<<0<<"\n";
				};
			};

			//rouResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<rouResidual[idx]<<"\n";			
			};

			//rovResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<rovResidual[idx]<<"\n";			
			};

			//rowResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<rowResidual[idx]<<"\n";			
			};

			//roResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<roResidual[idx]<<"\n";			
			};

			//roEResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<roEResidual[idx]<<"\n";			
			};

			//roU
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rou<<"\n";
			};

			//roV
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rov<<"\n";		
			};

			//dU_dx
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<gradCellsU[idx].x<<"\n";			
			};

			//dU_dy
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<gradCellsU[idx].y<<"\n";			
			};

			//Connectivity list for each cell			
			for (int i = 0; i<cells.size(); i++) {
				if (cells[i]->CGNSType == QUAD_4) {
				ofs<<toNaturalIndex[cells[i]->Nodes[0]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[1]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[2]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[3]]<<"\n";
				};
				if (cells[i]->CGNSType == TRI_3) {
				ofs<<toNaturalIndex[cells[i]->Nodes[0]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[1]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[2]]<<"\n";					
				};
			};
		};

		//3D case
		if (_grid.gridInfo.GridDimensions == 3) {

			//Access local nodes, faces, cells and flow data
			std::vector<Node*> nodes = _grid.nodes.getLocalNodes();
			std::vector<Face*> faces = _grid.faces.getLocalNodes();
			std::vector<Cell*> cells = _grid.cells.getLocalNodes();
			std::vector<ConservativeVariables*> data = U.getLocalNodes();

			ofs<<"VARIABLES= \"X\", \"Y\", \"Z\", \"Rho\", \"u\", \"v\", \"w\"";//, \"E\", \"T\", \"P\"";
			/*ofs<<", \"PStagnation\"";
			ofs<<", \"distance\"";
			ofs<<", \"rouResidual\"";
			ofs<<", \"rovResidual\"";
			ofs<<", \"rowResidual\"";
			ofs<<", \"roResidual\"";
			ofs<<", \"roEResidual\"";
			ofs<<", \"roU\"";
			ofs<<", \"roV\"";
			ofs<<", \"dU_dx\"";
			ofs<<", \"dU_dy\"";*/
			ofs<<"\n";
			
			ofs<<"ZONE T=\"D\"\n";
			ofs<<"N=" << _grid.nodes.size() << ", E=" << _grid.cells.size() <<", F=FEBLOCK, ";
			if (cells[0]->CGNSType == QUAD_4) {
				ofs<<"ET=QUADRILATERAL\n";
			};
			if (cells[0]->CGNSType == TRI_3) {
				ofs<<"ET=TRIANGLE\n";
			};
			if (cells[0]->CGNSType == HEXA_8) {
				ofs<<"ET=BRICK\n";
			};

			ofs<<"VARLOCATION = (NODAL, NODAL, NODAL, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED)\n";//, CELLCENTERED, CELLCENTERED, CELLCENTERED";
			/*ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<")\n";*/

			//Map all node global indexes to natural numbers
			std::map<int,int> toNaturalIndex;
			std::set<int> nodeIndexes = _grid.nodes.getAllIndexes();
			int counter = 1;
			for (std::set<int>::iterator it = nodeIndexes.begin(); it != nodeIndexes.end(); it++) toNaturalIndex[*it] = counter++;			

			//Nodes coordinates
			//X
			for (int i = 0; i<nodes.size(); i++) {
				ofs<<nodes[i]->P.x<<"\n";
			};

			//Y
			for (int i = 0; i<nodes.size(); i++) {
				ofs<<nodes[i]->P.y<<"\n";
			};

			//Z
			for (int i = 0; i<nodes.size(); i++) {
				ofs<<nodes[i]->P.z<<"\n";
			};

			////Solution data
			////Rho
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].ro<<"\n";
			};
				
			//u
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rou/U[idx].ro<<"\n";
			};
			//v
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rov/U[idx].ro<<"\n";		
			};
			//w
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].row/U[idx].ro<<"\n";
			};
			/*
			//E
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].roE/U[idx].ro<<"\n";
			};
			//T
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				double ro = U[idx].ro;
				double vx = U[idx].rou/U[idx].ro;
				double vy = U[idx].rov/U[idx].ro;
				double vz = U[idx].row/U[idx].ro;
				double E = U[idx].roE/U[idx].ro;
				double T = (E - (vx*vx+vy*vy+vz*vz)/2.0) / (medium.Cv);
				ofs<<T<<"\n";
			};
			//P
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				double ro = U[idx].ro;
				double vx = U[idx].rou/U[idx].ro;
				double vy = U[idx].rov/U[idx].ro;
				double vz = U[idx].row/U[idx].ro;
				double E = U[idx].roE/U[idx].ro;
				double P = (medium.Gamma-1.0) * ro * (E - (vx*vx+vy*vy+vz*vz)/2.0);
				ofs<<P<<"\n";
			};	

			//PStagnation
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				double ro = U[idx].ro;
				double vx = U[idx].rou/U[idx].ro;
				double vy = U[idx].rov/U[idx].ro;
				double vz = U[idx].row/U[idx].ro;
				double E = U[idx].roE/U[idx].ro;
				double PStagnation = (medium.Gamma-1.0) * U[idx].roE;
				ofs<<PStagnation<<"\n";
			};	

			//distance
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				if (_wallInfo.size() != 0) {
					ofs<<_wallInfo[idx].distance<<"\n";
				} else {
					ofs<<0<<"\n";
				};
			};

			//rouResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<rouResidual[idx]<<"\n";			
			};

			//rovResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<rovResidual[idx]<<"\n";			
			};

			//rowResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<rowResidual[idx]<<"\n";			
			};

			//roResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<roResidual[idx]<<"\n";			
			};

			//roEResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<roEResidual[idx]<<"\n";			
			};

			//roU
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rou<<"\n";
			};

			//roV
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rov<<"\n";		
			};

			//dU_dx
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<gradCellsU[idx].x<<"\n";			
			};

			//dU_dy
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<gradCellsU[idx].y<<"\n";			
			};*/

			//Connectivity list for each cell			
			for (int i = 0; i<cells.size(); i++) {
				if (cells[i]->CGNSType == QUAD_4) {
				ofs<<toNaturalIndex[cells[i]->Nodes[0]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[1]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[2]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[3]]<<"\n";
				};
				if (cells[i]->CGNSType == TRI_3) {
				ofs<<toNaturalIndex[cells[i]->Nodes[0]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[1]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[2]]<<"\n";					
				};
				if (cells[i]->CGNSType == HEXA_8) {
				ofs<<toNaturalIndex[cells[i]->Nodes[0]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[1]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[2]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[3]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[4]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[5]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[6]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[7]]<<"\n";					
				};
			};
		};

		ofs.close();
		return;

	};
	virtual void SaveToTechPlot(std::string fname) {
		std::ofstream ofs(fname);
		ofs<<std::scientific;
		if (_grid.gridInfo.GridDimensions == 1) {
			//Header
			ofs<<"VARIABLES= \"x\"\n";
			ofs<<"\"y\"\n";
			ofs<<"\"z\"\n";
			ofs<<"\"ro\"\n";
			ofs<<"\"u\"\n";
			ofs<<"\"v\"\n";
			ofs<<"\"w\"\n";
			ofs<<"\"P\"\n";
			ofs<<"\"E\"\n";
			ofs<<"\"T\"\n";
			ofs<<"ZONE\n";

			std::vector<Cell*> cells = _grid.cells.getLocalNodes();
			for (int i = 0; i<cells.size(); i++) {			
				Cell& c = *cells[i];								
				ofs<<c.CellCenter.x<<"\n";
				ofs<<c.CellCenter.y<<"\n";
				ofs<<c.CellCenter.z<<"\n";
				ofs<<U[c.GlobalIndex].ro<<"\n";
				ofs<<GetVelocity(U[c.GlobalIndex]).x<<"\n";
				ofs<<GetVelocity(U[c.GlobalIndex]).y<<"\n";
				ofs<<GetVelocity(U[c.GlobalIndex]).z<<"\n";
				ofs<<GetPressure(U[c.GlobalIndex])<<"\n";
				ofs<<U[c.GlobalIndex].roE/U[c.GlobalIndex].ro<<"\n";
				ofs<<GetTemperature(U[c.GlobalIndex])<<"\n";
				ofs<<"\n";
			};		
		};

		if (_grid.gridInfo.GridDimensions == 2) {
			//TO DO unify		
			//Header

			//Access local nodes, faces, cells and flow data
			std::vector<Node*> nodes = _grid.nodes.getLocalNodes();
			std::vector<Face*> faces = _grid.faces.getLocalNodes();
			std::vector<Cell*> cells = _grid.cells.getLocalNodes();
			std::vector<ConservativeVariables*> data = U.getLocalNodes();

			ofs<<"VARIABLES= \"X\", \"Y\", \"Rho\", \"u\", \"v\", \"w\", \"E\", \"T\", \"P\"";
			ofs<<", \"PStagnation\"";
			ofs<<", \"distance\"";
			ofs<<", \"rouResidual\"";
			ofs<<", \"rovResidual\"";
			ofs<<", \"rowResidual\"";
			ofs<<", \"roResidual\"";
			ofs<<", \"roEResidual\"";
			ofs<<", \"roU\"";
			ofs<<", \"roV\"";
			ofs<<"\n";
			
			ofs<<"ZONE T=\"D\"\n";
			ofs<<"N=" << _grid.nodes.size() << ", E=" << _grid.cells.size() <<", F=FEBLOCK, ";
			if (cells[0]->CGNSType == QUAD_4) {
				ofs<<"ET=QUADRILATERAL\n";
			};
			if (cells[0]->CGNSType == TRI_3) {
				ofs<<"ET=TRIANGLE\n";
			};

			ofs<<"VARLOCATION = (NODAL, NODAL, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<")\n";

			//Map all node global indexes to natural numbers
			std::map<int,int> toNaturalIndex;
			std::set<int> nodeIndexes = _grid.nodes.getAllIndexes();
			int counter = 1;
			for (std::set<int>::iterator it = nodeIndexes.begin(); it != nodeIndexes.end(); it++) toNaturalIndex[*it] = counter++;			

			//Nodes coordinates
			//X
			for (int i = 0; i<nodes.size(); i++) {
				ofs<<nodes[i]->P.x<<"\n";
			};

			//Y
			for (int i = 0; i<nodes.size(); i++) {
				ofs<<nodes[i]->P.y<<"\n";
			};

			////Solution data
			////Rho
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].ro<<"\n";
			};
				
			//u
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rou/U[idx].ro<<"\n";
			};
			//v
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rov/U[idx].ro<<"\n";		
			};
			//w
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].row/U[idx].ro<<"\n";
			};
			//E
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].roE/U[idx].ro<<"\n";
			};
			//T
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				double ro = U[idx].ro;
				double vx = U[idx].rou/U[idx].ro;
				double vy = U[idx].rov/U[idx].ro;
				double vz = U[idx].row/U[idx].ro;
				double E = U[idx].roE/U[idx].ro;
				double T = (E - (vx*vx+vy*vy+vz*vz)/2.0) / (medium.Cv);
				ofs<<T<<"\n";
			};
			//P
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				double ro = U[idx].ro;
				double vx = U[idx].rou/U[idx].ro;
				double vy = U[idx].rov/U[idx].ro;
				double vz = U[idx].row/U[idx].ro;
				double E = U[idx].roE/U[idx].ro;
				double P = (medium.Gamma-1.0) * ro * (E - (vx*vx+vy*vy+vz*vz)/2.0);
				ofs<<P<<"\n";
			};	

			//PStagnation
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				double ro = U[idx].ro;
				double vx = U[idx].rou/U[idx].ro;
				double vy = U[idx].rov/U[idx].ro;
				double vz = U[idx].row/U[idx].ro;
				double E = U[idx].roE/U[idx].ro;
				double PStagnation = (medium.Gamma-1.0) * U[idx].roE;
				ofs<<PStagnation<<"\n";
			};	

			//distance
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				if (_wallInfo.size() != 0) {
					ofs<<_wallInfo[idx].distance<<"\n";
				} else {
					ofs<<0<<"\n";
				};
			};

			//rouResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<rouResidual[idx]<<"\n";			
			};

			//rovResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<rovResidual[idx]<<"\n";			
			};

			//rowResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<rowResidual[idx]<<"\n";			
			};

			//roResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<roResidual[idx]<<"\n";			
			};

			//roEResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<roEResidual[idx]<<"\n";			
			};

			//roU
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rou<<"\n";
			};

			//roV
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rov<<"\n";		
			};

			//Connectivity list for each cell			
			for (int i = 0; i<cells.size(); i++) {
				if (cells[i]->CGNSType == QUAD_4) {
				ofs<<toNaturalIndex[cells[i]->Nodes[0]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[1]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[2]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[3]]<<"\n";
				};
				if (cells[i]->CGNSType == TRI_3) {
				ofs<<toNaturalIndex[cells[i]->Nodes[0]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[1]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[2]]<<"\n";					
				};
			};
		};

		//3D case
		if (_grid.gridInfo.GridDimensions == 3) {

			//Access local nodes, faces, cells and flow data
			std::vector<Node*> nodes = _grid.nodes.getLocalNodes();
			std::vector<Face*> faces = _grid.faces.getLocalNodes();
			std::vector<Cell*> cells = _grid.cells.getLocalNodes();
			std::vector<ConservativeVariables*> data = U.getLocalNodes();

			ofs<<"VARIABLES= \"X\", \"Y\", \"Z\", \"Rho\", \"u\", \"v\", \"w\"";//, \"E\", \"T\", \"P\"";
			/*ofs<<", \"PStagnation\"";
			ofs<<", \"distance\"";
			ofs<<", \"rouResidual\"";
			ofs<<", \"rovResidual\"";
			ofs<<", \"rowResidual\"";
			ofs<<", \"roResidual\"";
			ofs<<", \"roEResidual\"";
			ofs<<", \"roU\"";
			ofs<<", \"roV\"";
			ofs<<", \"dU_dx\"";
			ofs<<", \"dU_dy\"";*/
			ofs<<"\n";
			
			ofs<<"ZONE T=\"D\"\n";
			ofs<<"N=" << _grid.nodes.size() << ", E=" << _grid.cells.size() <<", F=FEBLOCK, ";
			if (cells[0]->CGNSType == QUAD_4) {
				ofs<<"ET=QUADRILATERAL\n";
			};
			if (cells[0]->CGNSType == TRI_3) {
				ofs<<"ET=TRIANGLE\n";
			};
			if (cells[0]->CGNSType == HEXA_8) {
				ofs<<"ET=BRICK\n";
			};

			ofs<<"VARLOCATION = (NODAL, NODAL, NODAL, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED)\n";//, CELLCENTERED, CELLCENTERED, CELLCENTERED";
			/*ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<")\n";*/

			//Map all node global indexes to natural numbers
			std::map<int,int> toNaturalIndex;
			std::set<int> nodeIndexes = _grid.nodes.getAllIndexes();
			int counter = 1;
			for (std::set<int>::iterator it = nodeIndexes.begin(); it != nodeIndexes.end(); it++) toNaturalIndex[*it] = counter++;			

			//Nodes coordinates
			//X
			for (int i = 0; i<nodes.size(); i++) {
				ofs<<nodes[i]->P.x<<"\n";
			};

			//Y
			for (int i = 0; i<nodes.size(); i++) {
				ofs<<nodes[i]->P.y<<"\n";
			};

			//Z
			for (int i = 0; i<nodes.size(); i++) {
				ofs<<nodes[i]->P.z<<"\n";
			};

			////Solution data
			////Rho
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].ro<<"\n";
			};
				
			//u
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rou/U[idx].ro<<"\n";
			};
			//v
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rov/U[idx].ro<<"\n";		
			};
			//w
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].row/U[idx].ro<<"\n";
			};
			/*
			//E
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].roE/U[idx].ro<<"\n";
			};
			//T
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				double ro = U[idx].ro;
				double vx = U[idx].rou/U[idx].ro;
				double vy = U[idx].rov/U[idx].ro;
				double vz = U[idx].row/U[idx].ro;
				double E = U[idx].roE/U[idx].ro;
				double T = (E - (vx*vx+vy*vy+vz*vz)/2.0) / (medium.Cv);
				ofs<<T<<"\n";
			};
			//P
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				double ro = U[idx].ro;
				double vx = U[idx].rou/U[idx].ro;
				double vy = U[idx].rov/U[idx].ro;
				double vz = U[idx].row/U[idx].ro;
				double E = U[idx].roE/U[idx].ro;
				double P = (medium.Gamma-1.0) * ro * (E - (vx*vx+vy*vy+vz*vz)/2.0);
				ofs<<P<<"\n";
			};	

			//PStagnation
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				double ro = U[idx].ro;
				double vx = U[idx].rou/U[idx].ro;
				double vy = U[idx].rov/U[idx].ro;
				double vz = U[idx].row/U[idx].ro;
				double E = U[idx].roE/U[idx].ro;
				double PStagnation = (medium.Gamma-1.0) * U[idx].roE;
				ofs<<PStagnation<<"\n";
			};	

			//distance
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				if (_wallInfo.size() != 0) {
					ofs<<_wallInfo[idx].distance<<"\n";
				} else {
					ofs<<0<<"\n";
				};
			};

			//rouResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<rouResidual[idx]<<"\n";			
			};

			//rovResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<rovResidual[idx]<<"\n";			
			};

			//rowResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<rowResidual[idx]<<"\n";			
			};

			//roResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<roResidual[idx]<<"\n";			
			};

			//roEResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<roEResidual[idx]<<"\n";			
			};

			//roU
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rou<<"\n";
			};

			//roV
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rov<<"\n";		
			};

			//dU_dx
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<gradCellsU[idx].x<<"\n";			
			};

			//dU_dy
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<gradCellsU[idx].y<<"\n";			
			};*/

			//Connectivity list for each cell			
			for (int i = 0; i<cells.size(); i++) {
				if (cells[i]->CGNSType == QUAD_4) {
				ofs<<toNaturalIndex[cells[i]->Nodes[0]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[1]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[2]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[3]]<<"\n";
				};
				if (cells[i]->CGNSType == TRI_3) {
				ofs<<toNaturalIndex[cells[i]->Nodes[0]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[1]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[2]]<<"\n";					
				};
				if (cells[i]->CGNSType == HEXA_8) {
				ofs<<toNaturalIndex[cells[i]->Nodes[0]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[1]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[2]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[3]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[4]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[5]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[6]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[7]]<<"\n";					
				};
			};
		};

		ofs.close();
		return;

	};
	virtual void SaveToTechPlotUndim(std::string fname, ConservativeVariables& var) {
		std::ofstream ofs(fname);
		ofs<<std::scientific;
		if (_grid.gridInfo.GridDimensions == 1) {
			//Header
			ofs<<"VARIABLES= \"x\"\n";
			ofs<<"\"ro\"\n";
			ofs<<"\"u\"\n";
			ofs<<"\"P\"\n";
			ofs<<"\"T\"\n";
			ofs<<"ZONE\n";

			std::vector<Cell*> cells = _grid.cells.getLocalNodes();
			for (int i = 0; i<cells.size(); i++) {			
				Cell& c = *cells[i];								
				ofs<<c.CellCenter.x<<"\n";
				ofs<<U[c.GlobalIndex].ro/var.ro<<"\n";
				ofs<<GetVelocity(U[c.GlobalIndex]).x/GetVelocity(var).x<<"\n";
				ofs<<GetPressure(U[c.GlobalIndex])/GetPressure(var)<<"\n";
				ofs<<GetTemperature(U[c.GlobalIndex])/GetTemperature(var)<<"\n";
				ofs<<"\n";
			};		
		};
		if (_grid.gridInfo.GridDimensions == 2) {
			//TO DO unify		
			//Header

			//Access local nodes, faces, cells and flow data
			std::vector<Node*> nodes = _grid.nodes.getLocalNodes();
			std::vector<Face*> faces = _grid.faces.getLocalNodes();
			std::vector<Cell*> cells = _grid.cells.getLocalNodes();
			std::vector<ConservativeVariables*> data = U.getLocalNodes();

			ofs<<"VARIABLES= \"X\", \"Y\", \"Rho\", \"u\", \"v\", \"w\", \"E\", \"T\", \"P\"";
			ofs<<", \"PStagnation\"";
			ofs<<", \"distance\"";
			ofs<<", \"rouResidual\"";
			ofs<<", \"rovResidual\"";
			ofs<<", \"rowResidual\"";
			ofs<<", \"roResidual\"";
			ofs<<", \"roEResidual\"";
			ofs<<", \"roU\"";
			ofs<<", \"roV\"";
			ofs<<", \"dU_dx\"";
			ofs<<", \"dU_dy\"";
			ofs<<"\n";
			
			ofs<<"ZONE T=\"D\"\n";
			ofs<<"N=" << _grid.nodes.size() << ", E=" << _grid.cells.size() <<", F=FEBLOCK, ";
			if (cells[0]->CGNSType == QUAD_4) {
				ofs<<"ET=QUADRILATERAL\n";
			};
			if (cells[0]->CGNSType == TRI_3) {
				ofs<<"ET=TRIANGLE\n";
			};

			ofs<<"VARLOCATION = (NODAL, NODAL, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<", CELLCENTERED";
			ofs<<")\n";

			//Map all node global indexes to natural numbers
			std::map<int,int> toNaturalIndex;
			std::set<int> nodeIndexes = _grid.nodes.getAllIndexes();
			int counter = 1;
			for (std::set<int>::iterator it = nodeIndexes.begin(); it != nodeIndexes.end(); it++) toNaturalIndex[*it] = counter++;			

			//Nodes coordinates
			//X
			for (int i = 0; i<nodes.size(); i++) {
				ofs<<nodes[i]->P.x<<"\n";
			};

			//Y
			for (int i = 0; i<nodes.size(); i++) {
				ofs<<nodes[i]->P.y<<"\n";
			};

			////Solution data
			////Rho
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].ro/var.ro<<"\n";
			};
				
			//u
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rou*var.ro/(U[idx].ro*var.rou)<<"\n";
			};
			//v
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rov*var.ro/(U[idx].ro*var.rov)<<"\n";		
			};
			//w
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].row*var.ro/(U[idx].ro*var.row)<<"\n";
			};
			//E
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].roE*var.ro/(U[idx].ro*var.roE)<<"\n";
			};
			//T
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				double ro = U[idx].ro/var.ro;
				double vx = U[idx].rou/(ro*var.rou);
				double vy = U[idx].rov/(ro*var.rov);
				double vz = U[idx].row/(ro*var.row);
				double E = U[idx].roE/(ro*var.roE);
				double T = (E - (vx*vx+vy*vy+vz*vz)/2.0) / (medium.Cv);
				ofs<<T<<"\n";
			};
			//P
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				double ro = U[idx].ro/var.ro;
				double vx = U[idx].rou/(ro*var.rou);
				double vy = U[idx].rov/(ro*var.rov);
				double vz = U[idx].row/(ro*var.row);
				double E = U[idx].roE/(ro*var.roE);
				double P = (medium.Gamma-1.0) * ro * (E - (vx*vx+vy*vy+vz*vz)/2.0);
				ofs<<P<<"\n";
			};	

			//PStagnation
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				double ro = U[idx].ro/var.ro;
				double vx = U[idx].rou/(ro*var.rou);
				double vy = U[idx].rov/(ro*var.rov);
				double vz = U[idx].row/(ro*var.row);
				double E = U[idx].roE/(ro*var.roE);
				double PStagnation = (medium.Gamma-1.0) * U[idx].roE;
				ofs<<PStagnation<<"\n";
			};	

			//distance
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				if (_wallInfo.size() != 0) {
					ofs<<_wallInfo[idx].distance<<"\n";
				} else {
					ofs<<0<<"\n";
				};
			};

			//rouResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<rouResidual[idx]<<"\n";			
			};

			//rovResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<rovResidual[idx]<<"\n";			
			};

			//rowResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<rowResidual[idx]<<"\n";			
			};

			//roResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<roResidual[idx]<<"\n";			
			};

			//roEResidual
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<roEResidual[idx]<<"\n";			
			};

			//roU
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rou<<"\n";
			};

			//roV
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<U[idx].rov<<"\n";		
			};

			//dU_dx
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<gradCellsU[idx].x<<"\n";			
			};

			//dU_dy
			for (int i = 0; i<cells.size(); i++) {
				int idx = cells[i]->GlobalIndex;
				ofs<<gradCellsU[idx].y<<"\n";			
			};

			//Connectivity list for each cell			
			for (int i = 0; i<cells.size(); i++) {
				if (cells[i]->CGNSType == QUAD_4) {
				ofs<<toNaturalIndex[cells[i]->Nodes[0]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[1]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[2]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[3]]<<"\n";
				};
				if (cells[i]->CGNSType == TRI_3) {
				ofs<<toNaturalIndex[cells[i]->Nodes[0]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[1]]<<" "
					<<toNaturalIndex[cells[i]->Nodes[2]]<<"\n";					
				};
			};
		};

		ofs.close();
		return;

	};
	void SaveCellGradientsToTechPlot2D(std::string fname) {
		std::ofstream ofs(fname);
		ofs<<std::scientific;

		//Access local nodes, faces, cells and flow data
		std::vector<Node*> nodes = _grid.nodes.getLocalNodes();
		std::vector<Face*> faces = _grid.faces.getLocalNodes();
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		std::vector<ConservativeVariables*> data = U.getLocalNodes();

		ofs<<"VARIABLES= \"X\", \"Y\", \"d_ro/dx\", \"d_ro/dy\", \"d_roU/dx\", \"d_roU/dy\", \"d_roV/dx\", \"d_roV/dy\"";
		ofs<<", \"d_roE/dx\"";
		ofs<<", \"d_roE/dy\"";
		ofs<<"\n";
			
		ofs<<"ZONE T=\"D\"\n";
		ofs<<"N=" << _grid.nodes.size() << ", E=" << _grid.cells.size() <<", F=FEBLOCK, ";
		if (cells[0]->CGNSType == QUAD_4) {
			ofs<<"ET=QUADRILATERAL\n";
		};
		if (cells[0]->CGNSType == TRI_3) {
			ofs<<"ET=TRIANGLE\n";
		};

		ofs<<"VARLOCATION = (NODAL, NODAL, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED";
		ofs<<", CELLCENTERED";
		ofs<<", CELLCENTERED";
		ofs<<")\n";

		//Map all node global indexes to natural numbers
		std::map<int,int> toNaturalIndex;
		std::set<int> nodeIndexes = _grid.nodes.getAllIndexes();
		int counter = 1;
		for (std::set<int>::iterator it = nodeIndexes.begin(); it != nodeIndexes.end(); it++) toNaturalIndex[*it] = counter++;			

		//Nodes coordinates
		//X
		for (int i = 0; i<nodes.size(); i++) {
			ofs<<nodes[i]->P.x<<"\n";
		};

		//Y
		for (int i = 0; i<nodes.size(); i++) {
			ofs<<nodes[i]->P.y<<"\n";
		};

		////Gradients data
		//d_ro/d_x
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<gradCellsRo[idx].x<<"\n";
		};
				
		//d_ro/d_y
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<gradCellsRo[idx].y<<"\n";
		};

		//d_roU/d_x
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<gradCellsRoU[idx].x<<"\n";	
		};

		//d_roU/d_y
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<gradCellsRoU[idx].y<<"\n";
		};

		//d_roV/d_x
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<gradCellsRoV[idx].x<<"\n";
		};

		//d_roV/d_y
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<gradCellsRoV[idx].y<<"\n";
		};

		//d_roE/d_x
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<gradCellsRoE[idx].x<<"\n";
		};	

		//d_roE/d_y
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<gradCellsRoE[idx].y<<"\n";
		};	

		//Connectivity list for each cell			
		for (int i = 0; i<cells.size(); i++) {
			if (cells[i]->CGNSType == QUAD_4) {
			ofs<<toNaturalIndex[cells[i]->Nodes[0]]<<" "
				<<toNaturalIndex[cells[i]->Nodes[1]]<<" "
				<<toNaturalIndex[cells[i]->Nodes[2]]<<" "
				<<toNaturalIndex[cells[i]->Nodes[3]]<<"\n";
			};
			if (cells[i]->CGNSType == TRI_3) {
			ofs<<toNaturalIndex[cells[i]->Nodes[0]]<<" "
				<<toNaturalIndex[cells[i]->Nodes[1]]<<" "
				<<toNaturalIndex[cells[i]->Nodes[2]]<<"\n";					
			};
		};
		ofs.close();
		return;
	};

	//Save solution slice in a file in Blazius coordinates
	void SaveSliceToTechPlot(std::string fname, double xStart, double Uinf, double xMin, double xMax, double yMin, double yMax) {
		std::ofstream ofs(fname);	
		ofs<<std::scientific;
		//Header
		ofs<<"VARIABLES= \"Y\"\n";
		ofs<<"\"U\"\n";
		ofs<<"\"Uright\"\n";
		ofs<<"ZONE\n";
		
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		std::vector<ConservativeVariables*> data = U.getLocalNodes();

		//Blazius data
		std::vector<double> yB;
		std::vector<double> uB;
		yB.push_back(0.0); uB.push_back(0.0);
		yB.push_back(0.1); uB.push_back(0.0664);
		yB.push_back(0.2); uB.push_back(0.1328);
		yB.push_back(0.3); uB.push_back(0.1989);
		yB.push_back(0.4); uB.push_back(0.2647);
		yB.push_back(0.5); uB.push_back(0.3298);
		yB.push_back(0.6); uB.push_back(0.3938);
		yB.push_back(0.7); uB.push_back(0.4563);
		yB.push_back(0.8); uB.push_back(0.5168);
		yB.push_back(0.9); uB.push_back(0.5748);
		yB.push_back(1.0); uB.push_back(0.6298);
		yB.push_back(1.1); uB.push_back(0.6813);
		yB.push_back(1.2); uB.push_back(0.7290);
		yB.push_back(1.3); uB.push_back(0.7725);
		yB.push_back(1.4); uB.push_back(0.8115);
		yB.push_back(1.5); uB.push_back(0.8460);
		yB.push_back(1.6); uB.push_back(0.8761);
		yB.push_back(1.7); uB.push_back(0.9018);
		yB.push_back(1.8); uB.push_back(0.9233);
		yB.push_back(1.9); uB.push_back(0.9411);
		yB.push_back(2.0); uB.push_back(0.9555);
		yB.push_back(2.1); uB.push_back(0.9670);
		yB.push_back(2.2); uB.push_back(0.9759);
		yB.push_back(2.3); uB.push_back(0.9827);
		yB.push_back(2.4); uB.push_back(0.9878);
		yB.push_back(2.5); uB.push_back(0.9915);
		yB.push_back(2.6); uB.push_back(0.9942);
		yB.push_back(2.7); uB.push_back(0.9962);
		yB.push_back(2.8); uB.push_back(0.9975);
		yB.push_back(2.9); uB.push_back(0.9984);		
		yB.push_back(5.0); uB.push_back(1.0);
		yB.push_back(std::numeric_limits<double>::max()); uB.push_back(1.0);
		for (int i = 0; i<yB.size(); i++) yB[i] *= 2;

		////Solution data
		////Rho
		for (int i = 0; i<cells.size(); i++) {			
			Cell& c = *cells[i];
			//Check if cell falls into zone
			if (c.CellCenter.x < xMin) continue;
			if (c.CellCenter.x > xMax) continue;
			if (c.CellCenter.y < yMin) continue;
			if (c.CellCenter.y > yMax) continue;

			//Compute blazius variables
			int idx = c.GlobalIndex;
			double y = c.CellCenter.y * sqrt(Uinf *U[idx].ro / ((c.CellCenter.x - xStart) * medium.Viscosity));
			double Unew = U[idx].rou / (U[idx].ro * Uinf);

			//interpolate blazius solution
			double Uright = 0;
			for (int j = 1; j<yB.size(); j++) {
				if ((y < yB[j]) && ( y >= yB[j-1] )) {
					Uright = uB[j-1] + (y - yB[j-1])*(uB[j] - uB[j-1])/(yB[j] - yB[j-1]);
					break;
				};
			};
			ofs<<y<<"\n";
			ofs<<Unew<<"\n";
			ofs<<Uright<<"\n";
			ofs<<"\n";
		};		

		ofs.close();
	};

	//Save to TecPlot boundary layer height
	void SaveBoundaryLayerHeightToTechPlot(std::string fname)
	{
		std::ofstream ofs(fname);
		for (std::map<int, FaceWallInfo>::iterator it = _wallFaces.begin(); it != _wallFaces.end(); it++)
		{
			int FaceInd = it->first;
			ofs << _grid.faces[FaceInd].FaceCenter.x << ' ' << _grid.faces[FaceInd].FaceCenter.y << ' ' << _grid.faces[FaceInd].FaceCenter.z << ' ';
			ofs << _boundaryLayerHeight[FaceInd] << '\n';
		};
		ofs.close();
	};

	//Load solution from TecPlot
	void LoadSolutionFromTecPlot2D(std::string fname) {
		std::ifstream ifs(fname);

		//variables for reading the file
		std::string str;
		double read;

		//start read
		std::getline(ifs, str);	//get first line with variables list
		std::getline(ifs, str); //get second line with Zone name
		std::getline(ifs, str);	//get third line with grid info
		std::getline(ifs, str); //get fours line with variables locations

		std::vector<Node*> nodes = _grid.nodes.getLocalNodes();
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();

		//read x and y coordinates of nodes
		for(int i=0; i<2*nodes.size(); i++) {
			ifs >> read;
		};

		//read density in all cells
		for(int i=0; i<cells.size(); i++) {
			ifs >> read;
			int idx = cells[i]->GlobalIndex;
			U[idx].ro = read;
		};
		//read u velocity component in all cells
		for(int i=0; i<cells.size(); i++) {
			ifs >> read;
			int idx = cells[i]->GlobalIndex;
			U[idx].rou = U[idx].ro*read;
		};
		//read v velocity component in all cells
		for(int i=0; i<cells.size(); i++) {
			ifs >> read;
			int idx = cells[i]->GlobalIndex;
			U[idx].rov = U[idx].ro*read;
		};
		//read w velocity component in all cells
		for(int i=0; i<cells.size(); i++) {
			ifs >> read;
			int idx = cells[i]->GlobalIndex;
			U[idx].row = U[idx].ro*read;
		};
		//read full energy in all cells
		for(int i=0; i<cells.size(); i++) {
			ifs >> read;
			int idx = cells[i]->GlobalIndex;
			U[idx].roE = U[idx].ro*read;
		};

		ifs.close();
	};

	//Save solution in a file
	void SaveSolution(std::string fname = "SavedSolution.sol")
	{
	 std::ofstream ofs;  
	 ofs.open(fname, std::ios::binary | std::ios::out);  
	 std::vector<Cell*> cells = _grid.cells.getLocalNodes();    
	 ofs.write( reinterpret_cast<char*>( &totalTime ), sizeof totalTime );
	 for(int i=0; i<cells.size(); i++)
		 {
		  int Index = cells[i]->GlobalIndex; //Global Index of cell
		  ofs.write( reinterpret_cast<char*>( &Index ), sizeof Index );
		  ofs.write( reinterpret_cast<char*>( &U[Index].ro ), sizeof U[Index].ro );
		  ofs.write( reinterpret_cast<char*>( &U[Index].roE ), sizeof U[Index].roE );
		  ofs.write( reinterpret_cast<char*>( &U[Index].rou ), sizeof U[Index].rou );
		  ofs.write( reinterpret_cast<char*>( &U[Index].rov ), sizeof U[Index].rov );
		  ofs.write( reinterpret_cast<char*>( &U[Index].row ), sizeof U[Index].row );   
		 };
	 ofs.close();
	};
	
	//Load a solution file
	void LoadSolution(std::string fname = "SavedSolution.txt")
	{
	 std::ifstream ifs;
	 ifs.open(fname, std::ios::binary | std::ios::in);
	 ifs.read( reinterpret_cast<char*>( &totalTime ), sizeof totalTime );
	 while(!ifs.eof())
	 {
		int Index;   
		ifs.read( reinterpret_cast<char*>( &Index ), sizeof Index );
		ifs.read( reinterpret_cast<char*>( &U[Index].ro ), sizeof U[Index].ro );
		ifs.read( reinterpret_cast<char*>( &U[Index].roE ), sizeof U[Index].roE );
		ifs.read( reinterpret_cast<char*>( &U[Index].rou ), sizeof U[Index].rou );
		ifs.read( reinterpret_cast<char*>( &U[Index].rov ), sizeof U[Index].rov );
		ifs.read( reinterpret_cast<char*>( &U[Index].row ), sizeof U[Index].row );   
	 };
	 ifs.close();
	};
	void ReadSolutionFromCGNS(std::string fname) {
		// TO DO process errors
		int fn;
		cg_open(fname.c_str(), CG_MODE_READ, &fn);

		// we know there is only one base (real working code would check!)
		int index_base=1;

		/* Check the cell and physical dimensions of the base. */	
		int physDim;
		int cellDim;
		char cgnsName[255];
		if(cg_base_read(fn, index_base, cgnsName, &cellDim, &physDim) != CG_OK) cg_error_exit();
		
		// we know there is only one zone (real working code would check!)
		int index_zone=1;
		// we know there is only one FlowSolution_t (real working code would check!)
		int index_flow=1;
		// get zone size (and name - although not needed here)
		cgsize_t sizes[3];
		char zone_name[255];
		cg_zone_read(fn, index_base, index_zone, zone_name, sizes);

		// check if size is consistent with grid size
		if (sizes[1] != U.size()) {
			std::cout <<"Error loading solution, wrong number of solution records!\n";
			return;
		};

		//   upper range index - use vertex dimensions
		//  checking GridLocation first (real working code would check
		//  to make sure there are no Rind cells also!):
		CGNS_ENUMT(GridLocation_t) location;
		char solution_name[255];
		cg_sol_info(fn,index_base,index_zone,index_flow, solution_name, &location);
		if (location != CellCenter) {
			std::cout <<"Error loading solution, GridLocation must be CellCenter!\n";
			return;
		};  
		// lower range index
		cgsize_t irmin[1];
		irmin[0] = 1;
		cgsize_t irmax[1];
		irmax[0] = U.size();
		// read flow solution
		// allocate input buffer for Double precision values
		double* inputBuffer = new double[U.size()];

		// load density
		cg_field_read(fn,index_base,index_zone,index_flow, "Density" , RealDouble, irmin, irmax, (void*)inputBuffer);
		for (int i = 0; i<U.size(); i++) {
			U[i+1].ro = inputBuffer[i];
		};

		// load velocity components
		cg_field_read(fn,index_base,index_zone,index_flow, "VelocityX" , RealDouble, irmin, irmax, (void*)inputBuffer);
		for (int i = 0; i<U.size(); i++) {
			U[i+1].rou = U[i+1].ro * inputBuffer[i];
		};

		// load Y component in case its not 1D problem
		for (int i = 0; i<U.size(); i++) U[i+1].rov = 0;			
		if (physDim > 1) {
			cg_field_read(fn,index_base,index_zone,index_flow, "VelocityY" , RealDouble, irmin, irmax, (void*)inputBuffer);
			for (int i = 0; i<U.size(); i++) {
				U[i+1].rov = U[i+1].ro * inputBuffer[i];
			};
		};

		// load Z component in case its not 1D or 2D problem
		for (int i = 0; i<U.size(); i++) U[i+1].row = 0;			
		if (physDim > 2) {
			cg_field_read(fn,index_base,index_zone,index_flow, "VelocityZ" , RealDouble, irmin, irmax, (void*)inputBuffer);
			for (int i = 0; i<U.size(); i++) {
				U[i+1].row = U[i+1].ro * inputBuffer[i];
			};
		};

		//load Velocities for axysymmetric problems
		if(IsAxySymmetric==true)
		{
			// load radial and axial velocity components
			cg_field_read(fn,index_base,index_zone,index_flow, "Axial_Velocity" , RealDouble, irmin, irmax, (void*)inputBuffer);
			for (int i = 0; i<U.size(); i++) {
				U[i+1].rou = U[i+1].ro * inputBuffer[i];
			};

			cg_field_read(fn,index_base,index_zone,index_flow, "Radial_Velocity" , RealDouble, irmin, irmax, (void*)inputBuffer);
			for (int i = 0; i<U.size(); i++) {
				U[i+1].rov = U[i+1].ro * inputBuffer[i];
			};

			//read tangential component
			cg_field_read(fn,index_base,index_zone,index_flow, "Swirl_Velocity" , RealDouble, irmin, irmax, (void*)inputBuffer);
			for (int i = 0; i<U.size(); i++) {
				U[i+1].row = U[i+1].ro * inputBuffer[i];
			};

		};

		// load stagnation energy
		cg_field_read(fn,index_base,index_zone,index_flow, "EnergyStagnation" , RealDouble, irmin, irmax, (void*)inputBuffer);  
		for (int i = 0; i<U.size(); i++) {
			U[i+1].roE = U[i+1].ro * inputBuffer[i];
		};

		// free buffer memory
		delete[] inputBuffer;

		// close CGNS file		
		cg_close(fn);
	};

	//New gradients section	
	//Compute scalar function gradient in each cell
	template<class T>
	void ComputeFunctionGradient( std::map<int, Vector>& grads, DistributedEntityManager<T>& values, double (Model<RiemannSolver>::*func)(const T&) ) {
		grads.clear();
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();

		//For each cell compute gradient of given function
		std::vector<Vector> nPoints;
		std::vector<double> nValues;
		for (int i = 0; i<cells.size(); i++) {			
			Cell& cell = *cells[i];			

			//Determine required set of point and values
			nPoints.clear();			
			nValues.clear();			

			//Add all neighbours
			for (int j = 0; j<cell.Faces.size(); j++) {
				Face& face = _grid.faces[cell.Faces[j]];				
				ConservativeVariables nU;
				if (face.isExternal) {
					nU = GetDummyCellValues(values[cell.GlobalIndex], face);
				} else {
					if (face.FaceCell_1 == cell.GlobalIndex) {
						nU = values[face.FaceCell_2];
					} else {
						nU = values[face.FaceCell_1];
					};
				};
				Vector nPoint;
				if (face.isExternal) {
					nPoint = 2 * (face.FaceCenter - cell.CellCenter) + cell.CellCenter;
				} else {
					if (face.FaceCell_1 == cell.GlobalIndex) {
						nPoint = _grid.cells[face.FaceCell_2].CellCenter;
					} else {
						nPoint = _grid.cells[face.FaceCell_1].CellCenter;
					};
				};
				double nValue = (this->*func)(nU);
				nPoints.push_back(nPoint);
				nValues.push_back(nValue);
			};

			double cellValue =  (this->*func)(U[cell.GlobalIndex]);
			grads[cell.GlobalIndex] = ComputeGradientByPoints(cell.CellCenter, cellValue, nPoints, nValues);
		};	

		return;
	};	

	//Computation of yPlus and uPlus dimensionless variables
	void ComputeWallVariables() {
		int sliceIndex;

		//TO DO Generalizi for curved surfaces
		//For each wall face compute
		for (std::pair<int, FaceWallInfo> pair : _wallFaces) {
			FaceWallInfo faceWI = pair.second;			
			if (faceWI.layerCells.size() == 0) {
				continue;
			};

			//TO DO Remove
			double xFace = _grid.faces[pair.first].FaceCenter.x;
			if ((xFace < 1.02) && (xFace > 1.0)) {
				sliceIndex = pair.first;
			};

			//Assume that first cell center falls into log region
			int firstCellIndex = faceWI.layerCells[0];
			const double E = 9.793;
			const double k = 0.4187;
			double y = _wallInfo[firstCellIndex].distance;
			double ro = GetDensity(U[firstCellIndex]);
			double u = GetVelocityX(U[firstCellIndex]);

			double A = log( E * y * ro / medium.Viscosity);
			double B = k * u;
			std::vector<double> params;
			params.push_back(A);
			params.push_back(B);

			//Compute shear velocity and shear stress
			alglib::real_1d_array x = "[1.0]";
			const double epsg = 1e-6;
			const double epsf = 0;
			const double epsx = 1e-6;
			alglib::ae_int_t maxits = 0;
			alglib::minlmstate state;
			alglib::minlmreport rep;
			
			alglib::minlmcreatev(1, x, 1e-6, state);
			alglib::minlmsetcond(state, epsg, epsf, epsx, maxits);
			alglib::minlmoptimize(state, LogWallRelation, NULL, &params);
			alglib::minlmresults(state, x, rep);

			//if (rep.terminationtype != 2) {
				//Error
				//throw Exception("");
			//};
			
			alglib::real_1d_array fCheck;
			fCheck.setlength(1);
			LogWallRelation(x, fCheck, &params);
			double checkF = fCheck[0];

			faceWI.shearVelocity = exp(x[0]);
			faceWI.shearStress = ro * faceWI.shearVelocity * faceWI.shearVelocity;							

			//In case yPlus < 10 assume laminar sublayer
			double yP = y * faceWI.shearVelocity * ro / medium.Viscosity;
			if ( yP < 20) {				
				faceWI.shearStress = medium.Viscosity * u / y;							
				faceWI.shearVelocity = sqrt(faceWI.shearStress / ro);
			};

			_wallFaces[pair.first] = faceWI;
		};

		//For each cell compute yPlus and uPlus
		for (std::pair<int, CellWallInfo> pair : _wallInfo) {
			CellWallInfo& cellWI = pair.second;
			double shearVelocity = _wallFaces[cellWI.wallFaceIndex].shearVelocity;
			double ro = GetDensity(U[pair.first]);
			//yPlus			
			cellWI.yPlus = cellWI.distance * shearVelocity * ro / medium.Viscosity;
			//uPlus
			cellWI.uPlus = GetVelocityX(U[pair.first]) / shearVelocity;
			_wallInfo[pair.first] = cellWI;
		};

		//TO DO Remove save profile
		//Output dimensionless profile for point around x = 1.0
		FaceWallInfo faceWI = _wallFaces[sliceIndex];
		std::vector<double> yPlus;
		std::vector<double> uPlus;
		function uDitrib(0);
		double Umax = 0;
		double RoAvg = 0;
		for (int cInd : faceWI.layerCells) {			
			RoAvg += GetDensity(U[cInd]);
			if (GetVelocityX(U[cInd]) > Umax) Umax = GetVelocityX(U[cInd]);
		};
		RoAvg /= faceWI.layerCells.size();
		uDitrib.value.push_back(0);
		uDitrib.x.push_back(0);
		for (int cInd : faceWI.layerCells) {
			double val = GetVelocityX(U[cInd])/Umax * (1.0 - GetVelocityX(U[cInd])/Umax);
			uDitrib.value.push_back(val);
			uDitrib.x.push_back(_grid.cells[cInd].CellCenter.y);
			uPlus.push_back(_wallInfo[cInd].uPlus);
			yPlus.push_back(_wallInfo[cInd].yPlus);
		};

		uDitrib.size = uDitrib.x.size();
		uDitrib.WriteData("111.dat");
		double deltaPP = uDitrib.TrapeziumIntegral();
		double RePP = Umax * deltaPP *RoAvg / medium.Viscosity;

		std::ofstream ofs("uPlus.dat");

		ofs<<std::scientific;
		//Header
		ofs<<"VARIABLES= \"yPlus\"\n";
		ofs<<"\"uPlus\"\n";		
		ofs<<"ZONE\n";

		for (int i = 0; i < uPlus.size(); i++) {
			ofs<<yPlus[i] << " " << uPlus[i] << "\n";
		};

		ofs.close();

		//Output skin friction coefficient along the plate
		ofs.open("CfGood.dat");

		ofs<<std::scientific;
		//Header
		ofs<<"VARIABLES= \"x\"\n";
		ofs<<"\"Cf\"\n";		
		ofs<<"\"CfTheory\"\n";		
		ofs<<"ZONE\n";

		for (std::pair<int, FaceWallInfo> pair : _wallFaces) {
			FaceWallInfo faceWI = pair.second;		
			int firstCellIndex = faceWI.layerCells[0];			
			double ro = GetDensity(U[firstCellIndex]);
			double Cf = faceWI.shearStress / ( 0.5 * Umax * Umax * RoAvg);
			double X = _grid.faces[pair.first].FaceCenter.x - 0.2;
			double ReX = X * Umax * RoAvg / medium.Viscosity;
			double CfT = 0.0592 * pow(ReX, -0.2);
			ofs<<X << " " << Cf << " " << CfT << "\n";
		};		

		ofs.close();

		function CfDistrib(0);		
		for (std::pair<int, FaceWallInfo> pair : _wallFaces) {
			FaceWallInfo faceWI = pair.second;		
			int firstCellIndex = faceWI.layerCells[0];			
			double ro = GetDensity(U[firstCellIndex]);
			double Cf = faceWI.shearStress / ( 0.5 * Umax * Umax * ro);
			CfDistrib.x.push_back(_grid.faces[pair.first].FaceCenter.x - 0.2);
			CfDistrib.value.push_back(Cf);
		};
		CfDistrib.size = CfDistrib.x.size();
		CfDistrib.WriteData("Cf.dat");
	};

};


#endif	